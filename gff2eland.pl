#!/usr/bin/env perl

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
) and (@ARGV or $ARGV{input}) or pod2usage( -verbose => 1 );

my ( $INH, $OUTH, $ERRH ) = _prepare_io( \%ARGV, \@ARGV );

while ( <$INH> ) {
    chomp;
    my @fields = split /\t/;

    my $mismatches = 'NM';
    my $target     = '';

    if ($fields[5] and $fields[7] ne q{.}) {
        my @mismatches = (0, 0, 0, 0); $mismatches[$fields[7]] = $fields[5];
        $mismatches = join ':', @mismatches;
        $target = "$fields[0]:$fields[3]" . ($fields[6] eq q{+} ? 'F' : 'R') . $fields[7];
    }
    my ($readid, $sequence) = $fields[2] =~ m{^(.+#\d/[12]):(.*)}g;

    die $fields[2] unless defined $readid and defined $sequence;

    print $OUTH join( "\t", $readid, $sequence, $mismatches, $target), "\n";
}


=head1 NAME

 gff2eland.pl - Convert GFF v3 to Eland v3

=head1 SYNOPSIS

 gff2eland.pl [OPTION]... [[-i] FILE]...

 gff2ekabd,pl your.gff -o your.eland3

=head1 DESCRIPTION

 Convert GFF v3 to Eland v3

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
