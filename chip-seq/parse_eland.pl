#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/sum/;

die "Usage: $0 <readsize> <lane> <eland_file>" unless @ARGV;

my $length  = shift @ARGV;
my $lane    = shift @ARGV;

while (<>) {
    my @eland_line = split /\t/, $_;

    next if $eland_line[10] =~ m/^QC$|NM$|^\d+:\d+:\d+$/;

    my ($seq_id, undef) = split /\./, $eland_line[10];
    my $coordinate = $eland_line[12];
    my $strand     = $eland_line[13] eq q{F} ? q{+} : q{-};
    my $pair_coord = $coordinate + $eland_line[19];

    my $center = 0;
    my $lib    = 0;
    if ($pair_coord >= $coordinate) {
        $center = ($pair_coord + $length - $coordinate) / 2 + $coordinate;
        $lib    =  $pair_coord + $length - $coordinate;
    }
    else {
        $center = ($coordinate + $length - $pair_coord) / 2 + $pair_coord;
        $lib    =  $coordinate + $length - $pair_coord;
    }

    next if $lib > 300;

    print join ("\t",
                $seq_id,
                q{.},
                "$lane",
                int $center,
                int $center,
                q{.},
                $strand,
                q{.},
                "/1=$coordinate;/2=$pair_coord;lib=$lib",
            ), "\n";
}
