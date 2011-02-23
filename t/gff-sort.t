#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use List::Util 'shuffle';

use Test::More qw(no_plan);

use FindBin;
use lib "$FindBin::Bin/../lib";
use GFF::Sort qw/gff_sort gff_is_sorted/;
use GFF::Parser;

my $orig = 't/sorted.gff';
my $dirty = 't/unsorted.gff';

for (1..10){
    shuffle_file($orig,$dirty);
    my $sorted = gff_sort(file => $dirty, overwrite => 0, column => 'start', manual => 0);
    ok(gff_is_sorted($sorted, 'start'), "file is sorted");
    ok(!gff_is_sorted($dirty, 'start'), "file is shuffled");
    unlink $sorted;
}

unlink $dirty;

sub shuffle_file{
    my ($in,$out) = @_;
    open my $inh, '<', $in;
    open my $outh, '>', $out;
    print $outh shuffle(<$inh>);
    close $outh;
    close $inh;
}
