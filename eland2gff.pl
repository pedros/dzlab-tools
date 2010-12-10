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
    'eland|e=s',
    _meta_options( \%ARGV ),
) and (@ARGV or $ARGV{input}) or pod2usage( -verbose => 1 );

my ( $INH, $OUTH, $ERRH ) = _prepare_io( \%ARGV, \@ARGV );


my $eland = index_eland( $ARGV{eland} );

print $ERRH scalar( keys %$eland), "\n";

my $n_reads;

while (<$INH>) {
    ++$n_reads and print $OUTH $_
    if exists $eland->{ get_id( $_, 2, -4, -3, -2 ) };
} 

print $ERRH $n_reads, "\n";

sub index_eland {
    my ($eland_file, %ids) = @_;

    open my $ELAND, '<', $eland_file or die "Can't open $eland_file: $!";

    while (my $eland_line = <$ELAND>) {
        my $id = get_id( $eland_line, 0, -3, -2, -1 );
        next unless $id;
        $ids{$id} = 1;
    }

    return \%ids;
}


sub get_id {
    my ($eland_line, $col, @parts) = @_;

    my $id = (split /\t/, $eland_line)[$col];
    $id    = join( ':', (split /:/, $id)[@parts] );

    return $id;
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
