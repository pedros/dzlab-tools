#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/max/;

my $cmp_file     = $ARGV[0] || 0;
my $width        = $ARGV[1] || 50;
my $lane         = $ARGV[2] || '';
my %col_windows  = ();
my %cmp_windows  = ();

while (<STDIN>) { # eg. cat 'file_to_window' | window.pl 'filtering_file' width
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
        $col_windows{$chr}{$i} = 0 if !exists $col_windows{$chr}{$i};
    }
}

if ($cmp_file) {
    open my $CMP_FILE, '<', $cmp_file or die "Can't open $cmp_file";
    while (<$CMP_FILE>) {
        my ($chr, $start, $end) = (split /\t/)[0,3,4];

        for (my $i = $start; $i < $end; $i += $width) {
            $cmp_windows{$chr}{$i} = [$start, $end];
        }
    }
}

for my $chr (sort {$a cmp $b} keys %col_windows) {

    for my $window (sort {$a <=> $b} keys %{$col_windows{$chr}}) {

        print join ("\t",
                    $chr,
                    q{.},
                    "${lane}_w$width",
                    $window,
                    $window + $width - 1,
                    $col_windows{$chr}{$window},
                    q{.},
                    q{.},
                    q{.},
                    #$cmp_windows{$chr}{$window}->[0],
                    #$cmp_windows{$chr}{$window}->[1],
                ), "\n" if !$cmp_file or exists $cmp_windows{$chr}{$window};
    }
}
