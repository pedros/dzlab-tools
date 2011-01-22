#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util;

my $gff_a;
my $gff_b;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'gff-a|a=s'  => \$gff_a,
    'gff-b|b=s'  => \$gff_b,
    'output|o=s' => \$output,
    'verbose|v'  => sub { use diagnostics; },
    'quiet|q'    => sub { no warnings; },
    'help|h'     => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'   => sub { pod2usage ( -verbose => 2 ); }
);

pod2usage ( -verbose => 1 )
unless $gff_a and $gff_b;

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

my $max = 0;
my %probes = ();

open my $GFF_A, '<', $gff_a or croak "Can't read $gff_a";
open my $GFF_B, '<', $gff_b or croak "Can't read $gff_b";
while (1) {
    my $tmp = <$GFF_A>;
    last if !defined $tmp;
    chomp $tmp;
    my $a = (split /\t/, $tmp)[5];

    my ($a_chr, $a_start, $a_end, $a_attr) = (split /\t/, $tmp)[0, 3, 4, 8];

    $tmp = <$GFF_B>;
    last if !defined $tmp;
    chomp $tmp;
    my $b = (split /\t/, $tmp)[5];

    my ($b_chr, $b_start, $b_end, $b_attr) = (split /\t/, $tmp)[0, 3, 4, 8];

    $probes{$a_chr}{$a_start}->[0] = $a_end;
    $probes{$b_chr}{$b_start}->[1] = $b_end;
    $probes{$a_chr}{$a_start}->[2] = $a;
    $probes{$b_chr}{$b_start}->[3] = $b;
    $probes{$a_chr}{$a_start}->[4] = $a_attr;
    $probes{$b_chr}{$b_start}->[5] = $b_attr;
}
close $GFF_A;
close $GFF_B;

for my $chr (sort keys %probes) {
    for my $start (sort {$a <=> $b} keys %{$probes{$chr}}) {

        my $a_end = $probes{$chr}{$start}->[0] || next;
        my $b_end = $probes{$chr}{$start}->[1] || next;

        croak "Non-matching windows: start:\t$start\tend_a:\t$a_end\tend_b:\t$b_end"
        if $a_end != $b_end;

        my $a    = $probes{$chr}{$start}->[2];
        my $b    = $probes{$chr}{$start}->[3];

        next if $a eq q{.} or $b eq q{.};

        my $frac = log_ratio ($a, $b);

        print join ("\t",
                    $chr,
                    q{.},
                    'a-b',
                    $start,
                    $a_end,
                    sprintf("%g", $frac),
                    q{.},
                    q{.},
                    "$probes{$chr}{$start}->[4];a=$a;b=$b\n",
                );
    }
}

sub log_ratio {
    my ($a, $b) = @_;

    croak "log_ratio requires two values."
    unless defined $a and defined $b;

    return $a - $b;
    return (log (1 / $b)) / log 2 if $a == 0 and $b != 0;
    return (log ($a / 1)) / log 2   if $a != 0 and $b == 0;
    return (log ($a / $b)) / log 2 if $a != 0 and $b != 0;
    return 0;
}




__END__


=head1 NAME

 normalize.pl -- Given two gff files with expression/instances # fields, computes log_2 ratio of each probe/window

=head1 SYNOPSIS

 ./normalize.pl --gff-a a-file.gff --gff-b b-file.gff > a-over-b.gff

=head1 DESCRIPTION

 Notes on computing log ratio: if window a or window b have value = 0, computes log ratio with value = 1 instead

=head1 OPTIONS

 -a, --gff-a   experiment GFF file
 -b, --gff-b   control GFF file
 -o, --output  output file

=head1 REVISION

 Version 0.0.1

 $Rev: 412 $:
 $Author: psilva $:
 $Date: 2010-10-08 11:28:40 -0700 (Fri, 08 Oct 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/normalize.pl $:
 $Id: normalize.pl 412 2010-10-08 18:28:40Z psilva $:

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
