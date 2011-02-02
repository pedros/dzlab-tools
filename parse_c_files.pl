#!/usr/bin/env perl
# ___UNDOCUMENTED___

use strict;
use warnings;

my $tissue = shift;
my $batch  = 'batch' . shift;
my $chr    = 'chr' . shift;


my $bp = 0;
if (-e "$ARGV[0].freq") {
    open my $FREQ, '<', "$ARGV[0].freq" or warn "Can't read $ARGV[0].freq";
    while (<$FREQ>) {
        if (m/'bp'/) {
            ($bp) = $_ =~ m/(\d+)/;
            last;
        }
    }
    close $FREQ;
}

my (%filtered, %unfiltered) = ();
while (<>) {
    chomp;

    next if ($_ =~ m/^#.*$|^\s*$/);

    my $filter = 0;

    my ($context, $attribute) = (split /\t/, $_)[2, -1];

    my ($c, $t)   = (split /;/, $attribute);

    $filter = 1 if $t =~ m/\*$/;

    ($c) = $c =~ m/(\d+)/;
    ($t) = $t =~ m/(\d+)/;

    $unfiltered{$context}{c} += $c;
    $unfiltered{$context}{t} += $t;

    $filtered{$context}{c} += $c unless $filter;
    $filtered{$context}{t} += $t unless $filter;
}

my $total_unfiltered_c = $unfiltered{CG}{c} + $unfiltered{CHG}{c} + $unfiltered{CHH}{c};
my $total_unfiltered_t = $unfiltered{CG}{t} + $unfiltered{CHG}{t} + $unfiltered{CHH}{t};

my $total_filtered_c = $filtered{CG}{c} + $filtered{CHG}{c} + $filtered{CHH}{c};
my $total_filtered_t = $filtered{CG}{t} + $filtered{CHG}{t} + $filtered{CHH}{t};

my $unfiltered_C_ratio = $total_unfiltered_c / ($total_unfiltered_c + $total_unfiltered_t);
my $filtered_C_ratio = $total_filtered_c / ($total_filtered_c + $total_filtered_t);

my $unfiltered_CG_ratio = $unfiltered{CG}{c} / (($unfiltered{CG}{c} + $unfiltered{CG}{t}));
my $filtered_CG_ratio   = $filtered{CG}{c}   / (($filtered{CG}{c}   + $filtered{CG}{t}));

my $unfiltered_CHG_ratio = $unfiltered{CHG}{c} / (($unfiltered{CHG}{c} + $unfiltered{CHG}{t}));
my $filtered_CHG_ratio   = $filtered{CHG}{c}   / (($filtered{CHG}{c}   + $filtered{CHG}{t}));

my $unfiltered_CHH_ratio = $unfiltered{CHH}{c} / (($unfiltered{CHH}{c} + $unfiltered{CHH}{t}));
my $filtered_CHH_ratio   = $filtered{CHH}{c}   / (($filtered{CHH}{c}   + $filtered{CHH}{t}));

print join ("\t",
            '#ID',
            'BP',
            'C',
            'CG',
            'CHG',
            'CHH',
            'T',
            'TG',
            'THG',
            'THH',
            'C_ratio',
            'CG_ratio',
            'CHG_ratio',
            'CHH_ratio',
            'filtered_C',
            'filtered_CG',
            'filtered_CHG',
            'filtered_CHH',
            'filtered_T',
            'filtered_TG',
            'filtered_THG',
            'filtered_THH',
            'filtered_C_ratio',
            'filtered_CG_ratio',
            'filtered_CHG_ratio',
            'filtered_CHH_ratio',
        ), "\n" if 0;

print join ("\t",
            join (q{|}, $tissue, $batch, $chr),
            $bp,
            $total_unfiltered_c,
            $unfiltered{CG}{c},
            $unfiltered{CHG}{c},
            $unfiltered{CHH}{c},
            $total_unfiltered_t,
            $unfiltered{CG}{t},
            $unfiltered{CHG}{t},
            $unfiltered{CHH}{t},
            $unfiltered_C_ratio,
            $unfiltered_CG_ratio,
            $unfiltered_CHG_ratio,
            $unfiltered_CHH_ratio,
            $total_filtered_c,
            $filtered{CG}{c},
            $filtered{CHG}{c},
            $filtered{CHH}{c},
            $total_filtered_t,
            $filtered{CG}{t},
            $filtered{CHG}{t},
            $filtered{CHH}{t},
            $filtered_C_ratio,
            $filtered_CG_ratio,
            $filtered_CHG_ratio,
            $filtered_CHH_ratio,
        ), "\n";
