#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

use Test::More qw(no_plan);
use FindBin;
use lib "$FindBin::Bin/../lib";
use Fasta qw/slurp_fasta bisulfite_convert/;

my $file = 't/test.fasta';
my $tmp = 't/test.fasta.tmp';

system("perl fasta_bsrc.pl $file > $tmp");

my $original = slurp_fasta($file);
my $converted = slurp_fasta($tmp);

#print Dumper $converted;

for my $seq (sort keys %$original) {
    my $original_forward = $original->{$seq};
    my $original_backward = $original_forward;
    $original_backward =~ tr/acgtACGT/tgcaTGCA/;

    $original_forward =~ tr/cC/tT/;
    $original_backward =~ tr/cC/tT/;
    $original_backward = reverse $original_backward;

    is($original_forward, $converted->{$seq});
    is($original_backward, $converted->{"RC_$seq"});
}

unlink $tmp;

