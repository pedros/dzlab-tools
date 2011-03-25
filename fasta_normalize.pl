#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';

use FindBin;
use lib "$FindBin::Bin/lib";
use Fasta;

my $file = shift // die "usage: $0 file.fasta > file.fasta.normalized";

my $fasta = slurp_fasta($file,normalize => 1);

for my $seq (sort keys %$fasta) {
    print format_fasta($seq, $fasta->{$seq});
}

