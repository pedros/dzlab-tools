#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

my $DATA_HANDLE = 'ARGV';
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

# the approximate value of the 97.5 percentile point
# of the normal distribution at 95% coverage
use constant NORMAL_SCORE => 1.96;

my @running_avg = ();
my @running_std = ();
my @running_var = ();
my @running_ste = ();
my @running_t95 = ();
my @num_scores  = ();

LOCUS:
while (<$DATA_HANDLE>) {

    chomp;

    my ($locus, @scores) = split /\t/;

  SCORE:
    for my $k (0 .. $scores - 1) {

        last SCORE if $k > @scores - 1;
        next SCORE if $scores[$k] =~ m/nan*/ixms;

        $running_avg[$k] = $scores[$k] unless defined $running_avg[$k];
        $running_std[$k] = $scores[$k] unless defined $running_ste[$k];
        $running_var[$k] = $scores[$k] unless defined $running_var[$k];
        $running_ste[$k] = $scores[$k] unless defined $running_ste[$k];
        $running_t95[$k] = $scores[$k] unless defined $running_t95[$k];
        $num_scores[$k]  = 0           unless defined $num_scores[$k];

        my $previous_avg = $running_avg[$k];

        $running_avg[$k]
        = $running_avg[$k]
        + ($scores[$k] - $running_avg[$k])
        / ++$num_scores[$k];

        $running_std[$k]
        = $running_std[$k]
        + ($scores[$k] - $previous_avg)
        * ($scores[$k] - $running_avg[$k]);

        $running_var[$k]
        = $running_std[$k]
        / ($num_scores[$k] - 1) if $num_scores[$k] > 1;

        $running_ste[$k]
        = sqrt ($running_var[$k])
        / sqrt ($num_scores[$k]);

        $running_t95[$k]
        = $running_ste[$k] * NORMAL_SCORE;
    }
}

print join ("\t", qw/bin mean mean-t95 mean+t95 std var ste scores/), "\n";

for my $k (0 .. $scores - 1) {

    print join ("\t",
                $k * $bin_width - int ($scores/2) * $bin_width,
                $running_avg[$k],
                $running_avg[$k] - $running_t95[$k],
                $running_avg[$k] + $running_t95[$k],
                sqrt ($running_var[$k]),
                $running_var[$k],
                $running_ste[$k],
                $num_scores[$k],
            ), "\n";

}

__END__
Running average, std and variance based on Knuth: TAOCP, Volume 2, p.~232
Mk = Mk-1+ (xk - Mk-1)/k
Sk = Sk-1 + (xk - Mk-1)*(xk - Mk).
Sk/(k - 1).

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
