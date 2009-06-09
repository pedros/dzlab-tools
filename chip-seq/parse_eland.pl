#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/sum/;

my $read_size;    # solexa sequences length
my $feature;      # third GFF field
my $pair_ends;    # for extracting (or not) necessary pair information
my $library_size; # expected value (for calculating center coordinates)
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'read-size|r=i'    => \$read_size,
    'feature|f=s'      => \$feature,
    'pair-ends|p'      => \$pair_ends,
    'library-size|l=i' => \$library_size,
    'output|o=s'       => \$output,
    'verbose|v'        => sub { use diagnostics; },
    'quiet|q'          => sub { no warnings; },
    'help|h'           => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'         => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $result and @ARGV and $read_size and $library_size;

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

while (<>) {

    # clean up and split input eland (export file)
    chomp;
    my @eland_line = split /\t/;
    next if
    $eland_line[10] =~ m/^QC$|NM$|^\d+:\d+:\d+$/;

    # parse eland line
    my ($seq_id, undef) = split /\./, $eland_line[10];
    my $coordinate = $eland_line[12];
    my $strand     = $eland_line[13] eq q{F} ? q{+} : q{-};

    my $center;    # center coordinate
    my $attribute; # Ninth GFF field

    if ($pair_ends) {
        my $pair_coord
        = $coordinate + $eland_line[19];  # see eland spec: field 19 holds the pair coordinate offset
        my $lib;                          # the observed library size, calculated from /1 start to /2 end

        if ($pair_coord >= $coordinate) { # equivalent to forward strand
            $center = ($pair_coord + $read_size - $coordinate) / 2 + $coordinate;
            $lib    =  $pair_coord + $read_size - $coordinate;
        }
        else {                            # equivalent to reverse strand
            $center = ($coordinate + $read_size - $pair_coord) / 2 + $pair_coord;
            $lib    =  $coordinate + $read_size - $pair_coord;
        }
        # minor quality control: don't allow observed library sizes larger than expected library sizes
        next if $lib > $library_size;

        $attribute = "/1=$coordinate;/2=$pair_coord;lib=$lib";
    }
    else { # if single ends, use read center coordinate
        $center = $coordinate + ($read_size / 2);
        $attribute = "/1=$coordinate"
    }

    # print new GFF line
    print join ("\t",
                $seq_id,
                q{.},
                "$feature",
                int $center,
                int $center,
                q{.},
                $strand,
                q{.},
                $attribute,
            ), "\n";
}


__END__


=head1 NAME

 parse_eland.pl - Convert Solexa export (eland) to GFF format

=head1 SYNOPSIS

 # convert eland alignment of 36bp reads (single-ends) with feature K27, library size 150 to gff
 # output to tmp, input file is s_5_export.txt (no switch)
 parse_eland.pl -r 36 -f K27-emb -l 150 -o tmp s_5_export.txt

=head1 DESCRIPTION

 Converts Solexa's eland format (aka s_*_export.txt) to GFF. Supports single and pair-ends files.

=head1 OPTIONS

 parse_eland.pl [OPTION]... [FILE]...

 -r, --read-size    Solexa sequences length (integer): REQUIRED
 -l, --library-size expected library size (for culling observed sizes beyond) REQUIRED
 -f, --feature      string to display in 3rd field in output GFF file
 -p, --pair-ends    input file contains information about pairs (pair coord. offset)
 -o, --output       filename to write results to (defaults to STDOUT)
 -v, --verbose      output perl's diagnostic and warning messages
 -q, --quiet        supress perl's diagnostic and warning messages
 -h, --help         print this information
 -m, --manual       print the plain old documentation page

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
