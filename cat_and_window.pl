#!/usr/bin/perl
package cat_and_window::Utils;    #all subs here
use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Spec;
use File::Path;
use File::Basename;

sub fasta_parse {        #why shouldn't this stop at underscores?
    my ($genome) = @_;
    my @groups;
    open my $GENOME, '<', $genome or croak "Can't open $genome: $!";
    while (<$GENOME>) { s/>(\S+)// and push @groups, $1; }
    tr/[A-Z]/[a-z]/ foreach @groups;
    close $GENOME or carp "Can't close $genome: $!";
    return @groups;
}

sub appender {
    my ( $file1, $file2 ) = @_;
    run_cmd("cat $file1 >> $file2");
}

sub window_gff_run {
    my ($file) = @_;
    run_cmd("window_gff.pl -f $file -w 1 -s 1 -o $file.merged");
    run_cmd("mv $file.merged $file");
}

sub window_gff_refactored_run {
    my ( $file1, $file2, $file3, $width, $step ) = @_;
    run_cmd(
"window_gff_REFACTORED.pl $file1 --width ${width} --step ${step} --absolute rice --no-skip --output $file2"
    );
    run_cmd("cat $file2 >> $file3");
}

sub run_cmd {
    my ($cmd) = @_;
    warn "-- CMD: $cmd\n";
    eval { system("$cmd") };
    croak "** failed to run command '$cmd': $@" if $@;
}


# something like this: "...batch1/single-c/Rice_Endorsperm-chr02.single-c-CHH.gff
#"${batch}/single[$sep]c/*[$sep]${group}[$sep]single[$sep]c*${context}*gff" 
my $files;

#make it more flexible. given a context, batch, root path, use the regex to find files. assume the same path format.
sub gimmefiles {    ##hmm...
    my ( $context, $batch, $regex ) = @_;
    $files = `find $batch -iname $context`;

}

sub build_regex_combinations {
    my @parts   = @_;
    my @regexes = ();

    ##TODO
    # http://stackoverflow.com/questions/127704/algorithm-to-return-all-combinations-of-k-elements-from-n

    return qr/(?:join q{|}, @regexes)/;
}


{
    my @matched_files = ();
    my $regex         = undef;

    sub get_files {
        my ( $root_dir, @parts ) = @_;

        return unless -d $root_dir;

        open my $ROOT_DIR, '<', $root_dir
        or croak "Can't open $root_dir: $!";

        use File::Spec;

        my @dir_contents
        = map { File::Spec->rel2abs($_) } readdir <$ROOT_DIR>;

        closedir $ROOT_DIR;

        $regex ||= build_regex_combinations (@parts);

        for (@dir_contents) {

            push @matched_files, $_
            if (File::Spec->splitpath($_))[2] =~ m/$regex/;

            my @child_files = get_files ($_, @parts)
            if -d $_;
        }

        return @matched_files;
    }
}


1;

package main;

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Spec;
use File::Path;
use File::Basename;

print STDERR "Starting cat and window run at ", `date`, "\n";

my $name;
my $genome;
my $contexts;
my @batches;
my $width    = 50;
my $step     = 50;
my $work_dir = File::Spec->curdir();
my $mode     = 'overwrite';            #other option is "append"
my $remove;
my $sep = q[.\s_:-];

my $result = GetOptions(
    'name|n=s'       => \$name,
    'genome|g=s'     => \$genome,
    'contexts|c=s'   => \$contexts,
    'batches|b=s{,}' => \@batches,
    'width|w=i'      => \$width,
    'step|s=i'       => \$step,
    'work_dir|d=s'   => \$work_dir,
    'mode|m=s'       => \$mode
);

my @groups = fasta_parse($genome);
my @contexts = split( " ", $contexts );


foreach my $batch (@batches) {$batch = "$work_dir/$batch";}

croak "Bad directory"
  unless -e "${work_dir}/post-processing/single-c"
      and -e "${work_dir}/post-processing/windows";

if ( $mode == 'overwrite' ) {
    run_cmd("rm -r ${work_dir}/post-processing/single-c");
    run_cmd("mkdir -p ${work_dir}/post-processing/single-c");
    run_cmd("rm -r ${work_dir}/post-processing/windows");
    run_cmd("mkdir -p ${work_dir}/post-processing/windows");
}

print STDERR "Concatenating single c files \n";

foreach my $group (@groups) {
    foreach my $batch (@batches) {
        foreach my $context (@contexts) {
            appender(
"${batch}/single[$sep]c/*[$sep]${group}[$sep]single[$sep]c*${context}*gff"                              
                ,    #this is where the regex stuff matters
"${work_dir}/post-processing/single-c/${name}_BS-Seq_${group}_${context}_w1_methylation.gff"
            );
        }
    }
    foreach my $context (@contexts) {
        window_gff_run(
"${work_dir}/post-processing/single-c/${name}_BS-Seq_${group}_${context}_w1_methylation.gff"
        );
    }
}
print STDERR "Done with code: $? \n";  # $? - what?

print STDERR "Merging and windowing concatenated single c files \n";

foreach my $group (@groups) {
    foreach my $context (@contexts) {
        window_gff_refactored_run(
"${work_dir}/post-processing/single-c/${name}_BS-Seq_${group}_${context}_w1_methylation.gff",
"${work_dir}/post-processing/windows/${name}_BS-Seq_${group}_${context}_w${width}_methylation.gff",
"${work_dir}/post-processing/windows/${name}_BS-Seq_all_${context}_w${width}_methylation.gff",
            $width,
            $step
        );
    }
}
print STDERR "Done with code: $? \n";

print STDERR "Finished cat and window run at ", `date`, "\n";

=head1 NAME

 cat_and window.pl - concatenate, merge, and window single c files

=head1 SYNOPSIS



=head1 DESCRIPTION

=head1 OPTIONS

 cat_and_window-pl [OPTION]... [FILE]...

 -m,   --mode
 -n,   --name
 -g,   --genome
 -c,   --contexts
 -b,   --batches   just give the name of each batch directory
 -w,   --width     (default 50)
 -s,   --step      (default 50)
 -d,   --work_dir  (default current directory)

=head1 REVISION

 

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
 Jon Graff <jonthegraff@berkeley.edu>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
