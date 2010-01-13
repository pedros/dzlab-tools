#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/reduce/;

my $output;
my $gff;

# Grabs and parses command line options
my $result = GetOptions (    
    'gff|g=s'     => \$gff,
    'output|o=s'  => \$output,
    'verbose|v'   => sub { use diagnostics; },
    'quiet|q'     => sub { no warnings; },
    'help|h'      => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'    => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

print "##gff-version 3\n";

while (<>) {
    chomp;
    my @fields    = split /\t/;
    my $attribute = "ID=$fields[0]; p=$fields[5]; maxwin=$fields[6]";

    if ($gff) {
        unless (ref $gff eq 'GLOB') {
            my $file = $gff;
            $gff     = undef;
            open $gff, '<', $file or croak "Can't open $file: $!"
        }

        my @gff_fields = ();
        while ( (split /\t/, <$gff>)[4] < $fields[2] ) {}
        until ( (my @gff_line = split /\t/, <$gff>)[3] > $fields[3] ) {
            push @gff_fields, \@gff_line;
        }

        my $max_window = reduce { $a->[5] > $b->[5] ? $a : $b } @gff_fields;
        $attribute .= "; maxstart=$max_window->[3]; maxend=$max_window->[4]";
    }

    print join ("\t",
                $fields[1],
                'chipotle',
                'locus',
                $fields[2],
                $fields[3],
                $fields[4],
                q{.},
                q{.},
                $attribute,
            ), "\n";
}


__END__


=head1 NAME

 chi2gff.pl -- Converts ChIPotle output files to gff format

=head1 SYNOPSIS

 # find maximum window coordinates in the original data, for each peak
 ./chi2gff.pl chipotle-out.peaks -o chipotle-out.gff -g original_data.gff

 -g is optional

=head1 DESCRIPTION


=head1 OPTIONS

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
