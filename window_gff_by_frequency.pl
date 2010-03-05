#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/max/;

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my $width = 50;
my $reference;
my $feature;
my $center;
my $no_skip;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'width|w=i'     => \$width,
    'feature|f=s'   => \$feature,
    'reference|r=s' => \$reference,
    'center|c'      => \$center,
    'no-skip|n'     => \$no_skip,
    'output|o=s'  => \$output,
    'verbose|v'   => sub { use diagnostics; },
    'quiet|q'     => sub { no warnings; },
    'help|h'      => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'    => sub { pod2usage ( -verbose => 2 ); }
);

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

my %col_windows = ();
while (<>) {
    my ($chr, $start, $end) = (split /\t/)[0,3,4];

    $start = $start + int (($end - $start + 1) / 2)
    if $center;

    if ($start % $width > 0) {
        $col_windows{$chr}{int ($start / $width) * $width + 1}++;
    }
    else {
        $col_windows{$chr}{$start + 1}++;
    }
}

$reference = index_fasta ($reference);

for my $chr (keys %col_windows) {

    my $last_coord;

    if (ref $reference eq 'HASH') {
        unless (exists $reference->{$chr}) {
            carp "Can't find $chr in fasta file. This chromosome will only have windows up to $last_coord";
            $last_coord = max keys %{$col_windows{$chr}};
        }
        else {$last_coord = $reference->{$chr}}
    }
    else {$last_coord = max keys %{$col_windows{$chr}};}
 
    if ($no_skip) {
        for (my $i = 1; $i <= $last_coord - $width + 1; $i += $width) {
            $col_windows{$chr}{$i} = 0 unless exists $col_windows{$chr}{$i};
        }
    }
}

for my $chr (sort {$a cmp $b} keys %col_windows) {

    for my $window (sort {$a <=> $b} keys %{$col_windows{$chr}}) {

        print join ("\t",
                    $chr,
                    q{.},
                    ${feature} // "w$width",
                    $window,
                    $window + $width - 1,
                    $col_windows{$chr}{$window},
                    q{.},
                    q{.},
                    q{.},
                ), "\n";
    }
}


sub index_fasta {
    my $reference_file = shift;

    my %reference = ();

    return unless $reference_file;

    # reads in the reference genome file into @fastaseq
    open my $REF, '<', "$reference_file" or croak "Can't open $reference for reading: $!";
    my @fastaseq = <$REF>;
    close $REF;

    # find and store indices for each chromosome change and corresponding descriptions
    my ( @idx, @dsc ) = ();
    for my $i ( 0 .. @fastaseq - 1 ) {
        if ( $fastaseq[$i] =~ m/^>/ ) {
            $fastaseq[$i] =~ s/>//g;
            $fastaseq[$i] = ( split /\s/, $fastaseq[$i] )[0];
            push @idx, $i;
            push @dsc, $fastaseq[$i];
        }
    }

    for my $j ( 0 .. @idx - 1 ) {
        my $line;
        if ( $j == scalar @idx - 1 ) {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. @fastaseq - 1]);
        }
        else {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. $idx[$j + 1] - 1]);
        }
        $line =~ s/[\n\r]//g;
        $reference{$dsc[$j]} = length $line;
    }
    return \%reference;
}


__END__


=head1 NAME

 gff_window_frequency.pl - Window GFF files by coordinate and measure frequency

=head1 SYNOPSIS

 # build windows until last coordinate in gff file
 gff_window_by_frequency.pl -w 50 -f 'some_feature' input_file.gff

 # build windows for whole genome
 gff_window_by_frequency.pl --width 50 --feature 'foo bar' --reference fasta.fa input_file.gff

=head1 DESCRIPTION

 Takes as input GFF files with single coordinate per line.
 Moves a sliding (non-overlapping) window across file.
 Measures frequency of coordinate occurrence per window.

=head1 OPTIONS

 gff_window_frequence.pl [OPTION]... [FILE]...

 -w, --width       sliding window width
 -f, --feature     third GFF feature
 -r, --reference   fasta file for computing chromosome lengths
 -c, --center      use center of regions to compute overlaps
 -n, --no-skip     output all windows, even those with no coverage
 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.2

 $Rev$:
 $Author$:
 $Date$:
 $HeadURL$:
 $Id$:

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
