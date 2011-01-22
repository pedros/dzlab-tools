#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/sum/;

my $DATA_HANDLE = 'ARGV';
my $output;
my $distance = 0;
my $score = 0;

# Grabs and parses command line options
my $result = GetOptions (
    'distance|d=f'=> \$distance,
    'score|s=f'   => \$score,
    'output|o=s'  => \$output,
    'verbose|v'   => sub { use diagnostics; },
    'quiet|q'     => sub { no warnings; },
    'help|h'      => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'    => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $result;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my @buffer = ();

while (<$DATA_HANDLE>) {
    next if /\s*#/;
    my @fields = split /\t/;

    if ( !@buffer || is_adjacent_by( \@buffer, \@fields, $distance, $score) ) {

        print join "\t", @fields[0..2], 1, $fields[3] - 1, 0, @fields[6..8]
        if !@buffer and $fields[3] > 1;
    
        push @buffer, \@fields;

    }

    else {
        flush_buffer( \@buffer );

        if ($buffer[-1][4] + 1 != $fields[4] and $fields[0] eq $buffer[-1][0]) {
            print join ("\t", @{$buffer[0]}[0..2], $buffer[-1][4] + 1, $fields[4] - 1, 0, @{$buffer[0]}[6..8])
        }

        @buffer = (\@fields);
    }
}

print join ("\t", @{$buffer[0]}[0..3], $buffer[-1]->[4], @{$buffer[0]}[5..8])
if @buffer;


sub is_adjacent_by {
    my ($buffer, $fields, $distance_range, $score_range) = @_;

    return
    (
        $fields->[0] eq $buffer->[-1][0]
        and abs( $fields->[3] - 1 - $buffer->[-1][4]) <= $distance_range
        and abs( $fields->[5] - 0 - $buffer->[-1][5]) <= $score_range
    );
}

sub flush_buffer {
    my ($buffer) = @_ or return;
    my @scores = map { $_->[5] } @$buffer;
    my $score = (sum @scores) / @scores;
    print join "\t", @{$buffer->[0]}[0..3], $buffer->[-1][4], sprintf( "%g", $score), @{$buffer->[0]}[6..8];
}

__END__


=head1 NAME

 collapse_gff.pl - Collapse GFF records to a single one if within a given distance and score range

=head1 SYNOPSIS

 # merge adjacent gff records with the same score
 perl collapse_gff.pl -o stuff_collapse.gff stuff_w1.gff

 # merge records separated by up to 10 bp and varying scores up to 2
 perl collapse_gff.pl -d 10 -s 2 -o stuff_collapse.gff stuff_w1.gff 

=head1 DESCRIPTION

=head1 OPTIONS

 collapse_gff.pl [OPTION]... [FILE]...

 -d, --distance-range max distance by which adjacent gff records can be separated for merging (default  0)
 -s, --score-range    max score difference by which adjacent gff records can be separated for merging (default  0)
 -o, --output         filename to write results to (defaults to STDOUT)
 -v, --verbose        output perl's diagnostic and warning messages
 -q, --quiet          supress perl's diagnostic and warning messages
 -h, --help           print this information
 -m, --manual         print the plain old documentation page

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
