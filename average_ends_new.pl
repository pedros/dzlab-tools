#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use Statistics::Descriptive;

my $output;
my $comments    = q{#};
my $scores      = 100;
my $bin_width   = 1;

# Grabs and parses command line options
my $result = GetOptions (
    'scores|s=i'    => \$scores,
    'bin-width|w=i' => \$bin_width,
    'output|o=s'    => \$output,
    'verbose|v'     => sub { use diagnostics; },
    'quiet|q'       => sub { no warnings; },
    'help|h'        => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'      => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result;

if ($output) {
    open my $USER_OUT, '>', $output
    or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my @stats = map { Statistics::Descriptive::Full->new() } (1 .. $scores);

LOCUS:
while (<ARGV>) {
    $_ =~ tr/\n\r//d;
    my ($locus, @scores) = split /\t/;

    for my $k (0 .. $scores - 1){
        next if $k >= @scores;
        my $s = $scores[$k];
        next if ($s eq 'na');
        $stats[$k]->add_data($s);
    }

}

print join ("\t", qw/bin mean std var ste numscores 25% 50% 75%/), "\n";

for my $k (0 .. $scores - 1){
    my $stat = $stats[$k];
    printf("%d\t" . ("%f\t" x 4) . "%d\t" . ("%f\t" x 3) . "\n",
        $k * $bin_width - int ($scores/2) * $bin_width,
        $stat->mean(),
        $stat->standard_deviation(),
        $stat->variance(),
        $stat->standard_deviation() / sqrt ($stat->count),
        $stat->count,
        quartiles($stat->get_data)
    );
}

sub quartiles{
    my @data = sort { $a<=>$b } @_;
    return @data[int(@data/4), int (@data/2), int (@data*3/4)];
}

=head1 NAME

 average_ends.pl - Average ends analysis bins wit statistics

=head1 SYNOPSIS

 # Providing --bin-width, -w lets the program output proper, absolute bin indices
 average_ends.pl --bin-width 100 --scores 100 input.ends -o output.dat

=head1 DESCRIPTION

 Expects an ends analysis file in which each record represents one locus, followed by a set number of scores per bin.
 A bin is a portion of the locus, for example from its 5' ends to its 3' end in the case of genes.

=head1 OPTIONS

 average_ends.pl [OPTION]... [FILE]...

 -s, --scores      Number of scores per locus record
 -w, --bin-width   Width of each scores bin for calculating absolute bin indices
 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: 249 $:
 $Author: psilva $:
 $Date: 2010-01-11 21:24:34 -0800 (Mon, 11 Jan 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/average_ends.pl $:
 $Id: average_ends.pl 249 2010-01-12 05:24:34Z psilva $:

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
