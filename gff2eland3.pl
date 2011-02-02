#!/usr/bin/env perl
# ___UNDOCUMENTED___

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use version; our $VERSION = qv('0.0.1');

use FindBin;
use lib "$FindBin::Bin/DZLab-Tools/lib";
use DZLab::Tools::RunUtils;

GetOptions(
    \%ARGV,
    'input|i=s', 'output|o=s', 'error|e=s',
    _meta_options( \%ARGV ),
) or pod2usage( -verbose => 1 );

my ( $INH, $OUTH, $ERRH ) = _prepare_io( \%ARGV, \@ARGV );

while ( <$INH> ) {
    

    #HWI-EAS105:2:1:8:808#0/1        +       chr5    10541462        GTTTTGATATTTAAAGTGAGTGGTGTNATATTAAGAAGAAGNATA   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII   0       26:A>N,41:G>N
    #chr5    U/U     HWI-EAS105:2:1:8:808#0/1:GTTTTGATATTTAAAGTGAGTGGTGTNATATTAAGAAGAAGNATA  10541462        10541506        1       +       2       target=CAGTCTCGATACCCAAAGTGAGTGGCGCAACATCAAGAAGAAGGATATC
    #chr3    U/NM    HWI-EAS172_628AB:4:1:1039:4511:CTTTATTCTTTTTTCTTATGTTGGTTTTGTCGTTTTTAAATGTTT  13719468        13719512        1       -       2       target=AAATCTATTCTTCTTCCTTATGTTGGTTTTGCCGTCTTCAAATGTACTC

    chomp;
    my @fields = split /\t/;
    my ($id, $sequence) = $fields[2] =~ m{(^.+:.+:.+:.+)?:(\w+)};
    my ($chr, $start, $strand, $mm, $alt) = @fields[0,3,6,7,5];
    next if $chr eq q{.};
    $mm = join( ',', (1) x $mm );

    print $OUTH join( "\t", $id, $strand, $chr, $start, $sequence, $sequence, $alt - 1, $mm ), "\n";
}

__DATA__


__END__

=head1 NAME

 APerlyName.pl - Short description

=head1 SYNOPSIS

 APerlyName.pl [OPTION]... [[-i] FILE]...

=head1 DESCRIPTION

 Long description

=head1 OPTIONS

 -i, --input       <string>     input filename                           (STDIN)
 -o, --output      <string>     output filename                          (STDOUT)
 -e, --error       <string>     output error filename                    (STDERR)
     --verbose     [integer]    print increasingly verbose error messages
     --quiet                    print no diagnostic or warning messages
     --version                  print current version
     --license                  print author's contact and copyright information
     --help                     print this information
     --manual                   print the plain old documentation page

=head1 VERSION

 0.0.1

=head1 REVISION

 $Rev: $:
 $Author: $:
 $Date: $:
 $HeadURL: $:
 $Id: $:

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
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
