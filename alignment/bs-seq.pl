#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Spec;
use File::Path;
use File::Basename;

# SIGINT trap. Ctrl-c triggers unclean exit.
$SIG{INT} = sub {croak "Received SIG$_[0]. Exiting...\n"};

my $left_read;   # required
my $right_read;  # required
my $reference;   # required
my $base_name    = 'out';
my $overwrite;
my $read_size    = 45;
my $var_length   = 0;
my $library_size = 300;
my $mismatches   = 2;
my $organism     = q{.};
my $batch        = 1;
my $window_size  = 50;
my $trust_dash_2 = 1;
my $single_ends  = 0;
my @left_splice;
my @right_splice;
my @groups       = ();

my @date = localtime (time);

my $out_dir = File::Spec->catdir (
    File::Spec->curdir(),
    sprintf "DZ_full-run_%4d-%02d-%02d_%02d.%02d.%02d",
    $date[5] + 1900, $date[4] + 1, $date[3], $date[2], $date[1], $date[0]
);

# Grabs and parses command line options
my $result = GetOptions (
    'left-read|l=s'        => \$left_read,
    'right-read|r=s'       => \$right_read,
    'reference|f=s'        => \$reference,
    'base-name|b=s'        => \$base_name,
    'overwrite|o'          => \$overwrite,
    'read-size|s=i'        => \$read_size,
    'variable-length|v=i'  => \$var_length,
    'library-size|k=i'     => \$library_size,
    'mismatches|n=i'       => \$mismatches,
    'organism|t=s'         => \$organism,
    'batch|i=i'            => \$batch,
    'window-size|w=i'      => \$window_size,
    'trust-dash-2|2'       => \$trust_dash_2,
    'single-ends|1=i'      => \$single_ends,
    'left-splice|ls=i{2}'  => \@left_splice,
    'right-splice|rs=i{2}' => \@right_splice,
    'groups|g=s{,}'        => \@groups,
    'out-directory|d=s'    => \$out_dir,
    'verbose|V'            => sub { use diagnostics; },
    'quiet|q'              => sub { no warnings; },
    'help|h'               => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'             => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $result and $left_read and $right_read and $reference;


# Get all chromosomes, pseudo-chromosomes, groups, etc, in the fasta reference file
# Discards all information after first blank character in fasta header
unless (@groups) {
    open my $REFERENCE, '<', $reference or croak "Can't open $reference: $!";

    while (<$REFERENCE>) {
        if (m/>/) {
            tr/A-Z/a-z/;
            m/>([^\s]+)/ && push @groups, $1;
        }
    }
    close $REFERENCE or carp "Can't close $reference: $!";
}

@left_splice  = (1, $read_size) unless @left_splice;
@right_splice = (1, $read_size) unless @right_splice;

# check whether output directory should be created, exists, should overwritten
if (-d $out_dir and $overwrite) {rmtree ( $out_dir, {keep_root => 1} )}
elsif (! -d $out_dir) {
    mkpath ( File::Spec->catfile ($out_dir, 'windows'),  {verbose => 1} );
    mkpath ( File::Spec->catfile ($out_dir, 'single-c'), {verbose => 1} );
}
#else {warn " overwrite $out_dir"}


my %files = (
    lfa   => File::Spec->catfile ($out_dir, File::Spec->canonpath($left_read))  . '.fa',
    lc2t  => File::Spec->catfile ($out_dir, File::Spec->canonpath($left_read))  . '.c2t',
    rfa   => File::Spec->catfile ($out_dir, File::Spec->canonpath($right_read)) . '.fa',
    rg2a  => File::Spec->catfile ($out_dir, File::Spec->canonpath($right_read)) . '.g2a',
    lel3  => File::Spec->catfile ($out_dir, File::Spec->canonpath($left_read))  . "_$left_splice[0]-$left_splice[1].eland3",
    rel3  => File::Spec->catfile ($out_dir, File::Spec->canonpath($right_read)) . "_$right_splice[0]-$right_splice[1].eland3",
    base  => File::Spec->catfile ($out_dir, $base_name)  . '.gff',
    log   => File::Spec->catfile ($out_dir, $base_name)  . '.log',
    split => _gen_files (File::Spec->catfile ($out_dir, $base_name), 'gff',  @groups),
    freq  => _gen_files (File::Spec->catfile ($out_dir, 'single-c', $base_name), 'single-c.gff', @groups),
    cont  => [_gen_files (File::Spec->catfile ($out_dir, 'single-c', $base_name), 'single-c-CG.gff', @groups),
              _gen_files (File::Spec->catfile ($out_dir, 'single-c', $base_name), 'single-c-CHG.gff', @groups),
              _gen_files (File::Spec->catfile ($out_dir, 'single-c', $base_name), 'single-c-CHH.gff', @groups)],
    wcont => [_gen_files (File::Spec->catfile ($out_dir, 'windows', $base_name), "w${window_size}-CG.gff", @groups),
              _gen_files (File::Spec->catfile ($out_dir, 'windows', $base_name), "w${window_size}-CHG.gff", @groups),
              _gen_files (File::Spec->catfile ($out_dir, 'windows', $base_name), "w${window_size}-CHH.gff", @groups)],
);

# convert reads
run_cmd ("fq_all2std.pl fq2fa $left_read > $files{lfa}")  unless file_exists($files{lfa});
run_cmd ("convert.pl c2t $files{lfa} > $files{lc2t}")     unless file_exists($files{lc2t});
unless ($single_ends) {
    run_cmd ("fq_all2std.pl fq2fa $right_read > $files{rfa}") unless file_exists($files{rfa});
    run_cmd ("convert.pl g2a $files{rfa} > $files{rg2a}")     unless file_exists($files{rg2a});
}

# convert genomes
run_cmd ("rcfas.pl $reference > $reference.rc")           unless file_exists("$reference.rc");
run_cmd ("convert.pl c2t $reference.rc > $reference.c2t") unless file_exists("$reference.c2t");
run_cmd ("convert.pl g2a $reference.rc > $reference.g2a") unless file_exists("$reference.g2a") or $single_ends;

# align with seqmap
run_cmd ("seqmap $mismatches $files{lc2t} $reference.c2t $files{lel3} /eland:3 /forward_strand /available_memory:8000 /cut:$left_splice[0],$left_splice[1]")   unless file_exists($files{lel3});
unless ($single_ends) {
    run_cmd ("seqmap $mismatches $files{rg2a} $reference.g2a $files{rel3} /eland:3 /forward_strand /available_memory:8000 /cut:$right_splice[0],$right_splice[1]") unless file_exists($files{rel3});
}
else {
    run_cmd ("seqmap $mismatches $files{lc2t} $reference.c2t $files{rel3} /eland:3 /forward_strand /available_memory:8000 /cut:$right_splice[0],$right_splice[1]") unless file_exists($files{rel3});
}

# get back original non-converted reads
run_cmd ("replace_reads.pl -f $files{lfa} -r $read_size --variable-length $var_length -s @left_splice  $files{lel3} > $files{lel3}.post") unless file_exists("$files{lel3}.post");
unless ($single_ends) {
    run_cmd ("replace_reads.pl -f $files{rfa} -r $read_size --variable-length $var_length -s @right_splice $files{rel3} > $files{rel3}.post") unless file_exists("$files{rel3}.post");
}
else {
    run_cmd ("replace_reads.pl -f $files{lfa} -r $read_size --variable-length $var_length -s @right_splice $files{rel3} > $files{rel3}.post") unless file_exists("$files{rel3}.post");
}

# make sure reads map together
run_cmd ("correlatePairedEnds.pl --left $files{lel3}.post --right $files{rel3}.post --reference $reference --output $files{base} --offset 0 --distance $library_size --readsize $read_size --trust-dash-2 $trust_dash_2 --single-ends $single_ends") unless file_exists($files{base});

# basic stats about the aligment
run_cmd ("collect_align_stats.pl $files{lel3}.post $files{rel3}.post $files{base} $organism $batch > $files{log}") unless file_exists($files{log});

# quantify methylation
for (@groups) {
    run_cmd ("split_gff.pl --sequence all $files{base}") unless (file_exists($files{split}->{$_}));
    run_cmd ("countMethylation.pl --ref $reference --gff $files{split}->{$_} --output $files{freq}->{$_} --sort") unless file_exists($files{freq}->{$_});
}

# window methylation counts into non-overlapping windows
for my $context (0..2) {
    for my $group (@groups) {
        run_cmd ("split_gff.pl --feature all $files{freq}->{$group}") unless file_exists($files{cont}->[$context]{$group});
        run_cmd ("window_gff.pl --gff-file $files{cont}->[$context]{$group} --width 1 --output $files{cont}->[$context]{$group}.merged") unless file_exists("$files{cont}->[$context]{$group}.merged");
        run_cmd ("window_gff.pl --gff-file $files{cont}->[$context]{$group}.merged --width $window_size --output $files{wcont}->[$context]{$group} --no-skip") unless file_exists($files{wcont}->[$context]{$group});
    }
}

### DONE


sub run_cmd {
    my ($cmd) = @_;
    warn "-- CMD: $cmd\n";
    eval {system ("$cmd")};
    croak "** failed to run command '$cmd': $@" if $@;
}

sub _gen_files {
    my ($base, $ext, @groups) = @_;

    my %split_files;
    for my $group (@groups) {
        $split_files{$group} = "$base-$group.$ext";
    }
    return \%split_files;
}

sub file_exists {
    my $file = shift;
    print "$file exists: skipping..." && return 1 if -f $file and -s $file;
    return 0;
}

__END__


=head1 NAME

 bs-seq.pl - Run dzlab's bisulfite sequencing analysis pipeline

=head1 SYNOPSIS

  bs-seq.pl --left-read s_7_1_sequence.txt --right-read s_7_2_sequence.txt --reference REF_COL.fa --base-name test-out --read-size 45 --library-size 300 --mismatches 2 --organism leco --batch 1 --left-splice 1 40 --right-splice 5 45

=head1 DESCRIPTION

=head1 OPTIONS

 bs-seq-pl [OPTION]... [FILE]...

 -l,  --left-read
 -r,  --right-read
 -f,  --reference
 -b,  --base-name
 -o,  --overwrite
 -s,  --read-size
 -v,  --variable-length
 -k,  --library-size
 -n,  --mismatches
 -t,  --organism
 -i,  --batch
 -ls, --left-splice
 -rs, --right-splice
 -1,  --single-ends
 -2,  --trust-dash-2
 -g,  --groups
 -d,  --out-directory
 -V,  --verbose     output perl's diagnostic and warning messages
 -q,  --quiet       supress perl's diagnostic and warning messages
 -h,  --help        print this information
 -m,  --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: $:
 $Author: $:
 $Date: $:
 $HeadURL: $:
 $Id: $:

=head1 AUTHOR

 Pedro Silva <psilva@nature.berkeley.edu/>
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
