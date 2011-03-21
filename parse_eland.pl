#!/usr/bin/env perl

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
my $eland_3;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'read-size|r=i'    => \$read_size,
    'feature|f=s'      => \$feature,
    'pair-ends|p'      => \$pair_ends,
    'library-size|l=i' => \$library_size,
    'eland-3|3'        => \$eland_3,
    'output|o=s'       => \$output,
    'verbose|v'        => sub { use diagnostics; },
    'quiet|q'          => sub { no warnings; },
    'help|h'           => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'         => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $result and @ARGV;

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

while (<>) {

    my $eland = ($eland_3 ?
              read_eland_3 ($_) :
              read_export ($_, $read_size, $library_size, $pair_ends)
          );

    print join ("\t",
                $eland->{chr},
                'el3',
                'read',
                $eland->{start},
                $eland->{end},
                $eland->{mm},
                $eland->{strand},
                q{.},
                q{.},
            ), "\n";

    # # print new GFF line
    # print join ("\t",
    #             $eland->{seq_id},
    #             q{.},
    #             ($feature ? $feature : 'parse_eland'),
    #             (exists $eland->{center} ? int $eland->{center} : $eland->{coordinate}),
    #             (exists $eland->{center} ? int $eland->{center} : $eland->{coordinate} + $read_size),
    #             q{.},
    #             $eland->{strand},
    #             q{.},
    #             (exists $eland->{attribute} ? $eland->{attribute} : "read=$eland->{read_id};seq=$eland->{sequence};mm=$eland->{mm}"),
    #         ), "\n";
} # done


sub read_eland_3 {
    my ($eland) = @_;
    chomp $eland;

    my ($id, $seq, $mm, $chr, $coord, $strand, $length) = split /\t/, $eland;

    carp $eland and return unless $chr and $seq and $mm and $chr;

    $chr =~ s/([^:]+)([RF])(\d)$//i;

    carp $eland and return unless defined $1 and defined $2 and defined $3;

    $coord  = $1;
    $strand = q{F} eq $2 ? q{+} : q{-};
    $mm     = $3;
    $chr   =~ s/://;

    return {
        id     => $id,
        start  => $coord,
        end    => $coord + length $seq,
        chr    => $chr,
        mm     => $mm,
        strand => $strand,
    };
}
# eland 3
# HWI-EAS105_0015:2:100:10008:17355#0/1   GTATGTGAATGTAAAGGATGTGGATGGTGTAGATGAATGTGTAGGAAGTGGATGGTGTAGATGACGAATGTCTAGG    0:0:1:0 chr3:873266F2
#                                                                                                                                 AT3G49470|chr3:18352908:18353474:246R1


sub read_export {
    my ($eland_line, $read_size, $library_size, $pair_ends)
    = @_;

    # clean up and split input eland (export file)
    chomp;
    my @eland_line = split /\t/;
    next if
    $eland_line[10] =~ m/^QC$|NM$|^\d+:\d+:\d+$/;

    # parse eland line

    my $read_id    = $eland_line[0];
    my $sequence   = $eland_line[8];
    my ($seq_id)   = split /\./, $eland_line[10], 1;
    my $coordinate = $eland_line[12];
    my $strand     = $eland_line[13] eq q{F} ? q{+} : q{-};
    my $mismatch   = 0;
    #(defined $eland_line[14] ? $eland_line[14] =~ tr/A-Za-z// : 0;

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

    return {
        read_id    => $read_id,
        sequence   => $sequence,
        mm         => $mismatch,
        strand     => $strand,
        coord      => $coordinate,
        center     => $center,
        attribute  => $attribute,
    }
}



__END__


=head1 NAME

 parse_eland.pl - Convert Solexa export (eland) to GFF format. strips off read ID's, extrapolates . Each line is a read. 
 Score is number of mismatches a read has to the region.

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

 $Rev: 377 $:
 $Author: psilva $:
 $Date: 2010-07-08 13:50:33 -0700 (Thu, 08 Jul 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/parse_eland.pl $:
 $Id: parse_eland.pl 377 2010-07-08 20:50:33Z psilva $:

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
