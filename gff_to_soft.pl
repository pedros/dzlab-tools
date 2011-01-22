#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

# Check required command line parameters
unless (@ARGV > 0) {
    pod2usage ( -verbose => 1 );
}

# Grabs and parses command line options
my $result = GetOptions (
    'verbose|v' => sub { use diagnostics; },
    'quiet|q'   => sub { no warnings; },
    'help|h'    => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'  => sub { pod2usage ( -verbose => 2 ); }
);

my $ID_REF = 0;

print "#ID_REF = \n";
print "#VALUE = Fractional methylation: c / (c + t)\n";
print "#NoC = Number of sequenced methylated cytosines\n";
print "#NoT = Number of sequenced non-methylated cytosines\n";
print "Coord = Forward strand position\n";
print "Strand = Orientation of alignment\n";
print "!Sample_table_begin\n";
print join "\t", qw(ID_REF VALUE NoC NoT Chr Coord Strand), "\n";

while (<>) {
    next if ($_ =~ m/^#.*$|^\s*$/);
    my @fields = split /\t/, $_;
    my ($c, $t) = split(/;/, $fields[-1]);
    ($c) = $c =~ m/(\d+)/;
    ($t) = $t =~ m/(\d+)/;

    print join ("\t",
                ++$ID_REF,
                sprintf("%g", $fields[5]),
                $c,
                $t,
                $fields[0],
                $fields[3],
                $fields[6],
                "\n"
            );
}

print '!Sample_table_end';

__END__


=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

=head1 REVISION

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
