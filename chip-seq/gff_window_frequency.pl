#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/max/;

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my $width = 50;
my $feature;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'width|w=i'   => \$width,
    'feature|f=s' => \$feature,
    'output|o=s'  => \$output,
    'verbose|v'   => sub { use diagnostics; },
    'quiet|q'     => sub { no warnings; },
    'help|h'      => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'    => sub { pod2usage ( -verbose => 2 ); }
);

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

my %col_windows = ();
while (<>) {
    my ($chr, $start) = (split /\t/)[0,3];

    if ($start % $width > 0) {
        $col_windows{$chr}{int ($start / $width) * $width + 1}++;
    }
    else {
        $col_windows{$chr}{$start + 1}++;
    }
}

for my $chr (keys %col_windows) {
    my $last_coord = max keys %{$col_windows{$chr}};

    for (my $i = 1; $i < $last_coord; $i += $width) {
        $col_windows{$chr}{$i} = 0 unless exists $col_windows{$chr}{$i};
    }
}


for my $chr (sort {$a cmp $b} keys %col_windows) {

    for my $window (sort {$a <=> $b} keys %{$col_windows{$chr}}) {

        print join ("\t",
                    $chr,
                    q{.},
                    "${feature}_w$width",
                    $window,
                    $window + $width - 1,
                    $col_windows{$chr}{$window},
                    q{.},
                    q{.},
                    q{.},
                ), "\n";
    }
}


__END__


=head1 NAME

 gff_window_frequency.pl - Window GFF files by coordinate and measure frequency

=head1 SYNOPSIS

 gff_window_frequency.pl -w 50 -f 'some_feature' input_file.gff

=head1 DESCRIPTION

 Takes as input GFF files with single coordinate per line.
 Moves a sliding (non-overlapping) window across file.
 Measures frequency of coordinate occurrence per window.

=head1 OPTIONS

 gff_window_frequence.pl [OPTION]... [FILE]...

 -w, --width       sliding window width
 -f, --feature     third GFF file
 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

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
