#!/usr/bin/env perl

package cat_and_window;

use Exporter;
BEGIN {
    our @ISA = qw(Exporter);
    our @EXPORT = qw(fasta_parse appender window_gff_run window_gff_refactored_run run_cmd build_regex_combinations get_files get_file os_indie);
}

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Spec;
use File::Path;
use File::Basename;
use constant TEST => 0;   #modify the necessary subs so that if TEST is 1, they only print out what they're doing, but don't actually do it.

sub fasta_parse {
    my ($genome) = @_;
    my @groups;
    open my $GENOME, '<', $genome or croak "Can't open $genome: $!";
    while (<$GENOME>) { s/>(\S+)// and push @groups, lc $1; }
    close $GENOME or carp "Can't close $genome: $!";
    return @groups;
}

sub appender {
    my ( $file1, $file2 ) = @_;    
    if ( TEST ) {
	print STDERR "Appending $file1 to $file2\n";
	open my $FILE2, '>', $file2 or croak "Can't open $file2: $!";
	close $FILE2;
    }
    else {
	open my $FILE1, '<',  $file1 or croak "Can't open $file1: $!";
	open my $FILE2, '>>', $file2 or croak "Can't open $file2: $!";
	while (<$FILE1>) { print $FILE2 "$_"; }
	close $FILE1 or croak "Can't close $file1: $!";
	close $FILE2 or croak "Can't close $file2: $!";
    }
}

sub window_gff_run {
    my ($file) = @_;
    if ( TEST ) {
	print STDERR "Running window_gff on $file\n"
    }
    else {
    run_cmd("window_gff.pl -f $file -w 1 -s 1 -o $file.merged");
    rename "$file.merged", $file;
    }
}

sub window_gff_refactored_run {
    my ( $file1, $file2, $file3, $width, $step ) = @_;
    if ( TEST ) {
	print STDERR "Running window_gff_refactored on $file1 with width of $width and step of $step. 
Outputting to $file2, then appending $file2 to $file3\n\n";
    }
    else {
	run_cmd(
	    "window_gff_REFACTORED.pl $file1 --width ${width} --step ${step} --absolute rice --no-skip --output $file2"
	    );
	appender( $file2, $file3 );
    }
}

sub run_cmd {
    my ($cmd) = @_;
    warn "-- CMD: $cmd\n";
    eval { system("$cmd") };
    croak "** failed to run command '$cmd': $@" if $@;
}

sub build_regex_combinations {
    my @parts   = @_;
    my @regexes = ();
    use Algorithm::Permute qw(permute);
    Algorithm::Permute::permute { push @regexes, join '.*', @parts; } @parts;
    my $regex = join q{|}, @regexes;
    return qr/$regex/i;
}

{
    my @matched_files = ();
    my %regexes       = ();
    my %matched_files = ();

    sub get_files {
        my ( $root_dir, @parts ) = @_;

        return unless -d $root_dir;

        my $normal_params = join q{:}, $root_dir, @parts;
        return @{ $matched_files{$normal_params} }
          if exists $matched_files{$normal_params};

        opendir my $ROOT_DIR, $root_dir
	    or croak "Can't open $root_dir: $!";

 #        my @dir_contents = map { File::Spec->rel2abs($_) }         rel2abs doesn't work when you're in another directory.
#	grep { $_ !~ m/^\.+$/ } readdir $ROOT_DIR;      

	my @dir_contents = grep { $_ !~ m/^\.+$/ } readdir $ROOT_DIR;
	foreach (@dir_contents) {
	    $_ = "$root_dir/$_";}

	closedir $ROOT_DIR;

        my $regex = $regexes{ join q{:}, @parts } //=
          build_regex_combinations(@parts);

        for (@dir_contents) {

            push @{ $matched_files{$normal_params} }, $_
              if ( File::Spec->splitpath($_) )[2] =~ m/$regex/;

            my @child_files = get_files( $_, @parts )
               if -d $_;
        }
      
        return 'ARRAY' eq ref $matched_files{$normal_params}
	? @{ $matched_files{$normal_params} }
	: ();
    }
}

sub get_file
{ #Calls get_files. Finds the single file matching @parts in $root_dir. if there are more than one files, croaks. If there are no files, returns $default_file.
    my ( $root_dir, $default_file, @parts ) = @_
      ; #To indicate that there must be exactly one file, input $default_file as "no default".
    my @files = get_files( $root_dir, @parts );
    if ( @files > 1 ) {
        croak
	    "Expected only one file matching @parts in $root_dir . Instead there are ",
	    scalar @files, " that match.";
    }
    if ( 0 == @files ) {
        if ( $default_file =~ m/no.?default/i ) {
            croak
"Expected exactly one file matching @parts in $root_dir . Instead there are ",
              scalar @files, " that match.";
        }
        $files[0] = $default_file;
    }
    return $files[0];
}

sub os_indie {
    my ($filepath) = @_;
    my @filepath = split m'/', $filepath;
    return File::Spec->catfile (@filepath);
}

1;

package main;

import cat_and_window;
use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Spec;
use File::Path;
use File::Basename;

my $name;
my $genome;
my @contexts;
my @batches;    
my $width    = 50;
my $step     = 50;
my $work_dir = File::Spec->curdir();
my $append = 1;
my $overwrite = 0;

my $result = GetOptions(
    'name|n=s'       => \$name,
    'genome|g=s'     => \$genome,
    'contexts|c=s{,}'   => \@contexts,
    'batches|b=s{,}' => \@batches,
    'width|w=i'      => \$width,
    'step|s=i'       => \$step,
    'work_dir|d=s'   => \$work_dir,
    'overwrite|o'    => sub {$overwrite = 1; $append = 0;}
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $name and $genome and @contexts and @batches and ($overwrite xor $append);

my @groups = fasta_parse($genome);
#my @contexts = split( " ", $contexts );
$work_dir = File::Spec->rel2abs($work_dir);

foreach my $batch (@batches) { $batch = os_indie ("$work_dir/$batch"); }

#croak "Bad directory"
 # unless -d os_indie("${work_dir}/post-processing/single-c")
  #    and -d os_indie("${work_dir}/post-processing/windows");

print STDERR "Starting cat and window run at ", `date`, "\n";

if ( $overwrite ) {
    rmtree os_indie("${work_dir}/post-processing/single-c");
    rmtree os_indie("${work_dir}/post-processing/windows");
}

mkpath os_indie("${work_dir}/post-processing/single-c");
mkpath os_indie("${work_dir}/post-processing/windows");

print STDERR "Concatenating single c files \n";

foreach my $group (@groups) {
    foreach my $batch (@batches) {
        foreach my $context (@contexts) {
            my $sc_file =
              get_file( os_indie ("$batch/single-c"), "no default", $group, $context,
                "single-c" );
            my $pp_file = get_file(
                os_indie ("{$work_dir}/post-processing/single-c"),
		os_indie ("${work_dir}/post-processing/single-c/${name}_BS-Seq_${group}_${context}_w1_methylation.gff"),
                $group,
                $context,
                $name
            );
	    
            appender( $sc_file, $pp_file );
        }
    }

    foreach my $context (@contexts) {
        my $pp_file = get_file( os_indie ("$work_dir/post-processing/single-c"),
            "no default", $group, $context, $name );
        window_gff_run($pp_file);
    }
}
print STDERR "Done with code: $? \n";

print STDERR "Merging and windowing concatenated single c files \n";

foreach my $group (@groups) {
    foreach my $context (@contexts) {
        my $pp_sc_file = get_file( os_indie ("$work_dir/post-processing/single-c"),
            "no default", $group, $context, $name );
        my $pp_w_file = get_file(
            os_indie ("$work_dir/post-processing/windows"),
	    os_indie ("${work_dir}/post-processing/windows/${name}_BS-Seq_${group}_${context}_w${width}_methylation.gff"),
            $group,
            $context,
            $name,
            $width
        );
        my $pp_w_all_file = get_file(
            os_indie ("$work_dir/post-processing/windows"),
	    os_indie ("${work_dir}/post-processing/windows/${name}_BS-Seq_all_${context}_w${width}_methylation.gff"),
            "all",
            $context,
            $name,
            $width
        );

        window_gff_refactored_run( $pp_sc_file, $pp_w_file, $pp_w_all_file,
            $width, $step );
    }
}
print STDERR "Done with code: $? \n";

print STDERR "Finished cat and window run at ", `date`, "\n";

=head1 NAME

 cat_and window.pl - concatenate, merge, and window single c files

=head1 SYNOPSIS

usage:
perl /home/jgraff/workspace/bisulfite/trunk/cat_and_window.pl -n Arabidopsis_BSWT -g /work/genomes/AT/TAIR_reference.fas -c CG CHG CHH -b batch-1 batch-2 -d /work/Arabidopsis_BS_devin/WT


=head1 DESCRIPTION

=head1 OPTIONS

 cat_and_window-pl [OPTION]... [FILE]...

 -n,   --name
 -g,   --genome
 -c,   --contexts
 -b,   --batches   just give the name of each batch directory, like "batch-1"
 -w,   --width     (default 50)
 -s,   --step      (default 50)
 -d,   --work_dir  (default current directory) 
 -o,   --overwrite    The default mode is append. This changes it to overwrite.


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
