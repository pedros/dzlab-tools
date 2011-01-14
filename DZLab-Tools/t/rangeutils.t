#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use List::Util qw/shuffle/;

use DZLab::Tools::RangeUtils;

#use Test::Simple qw(no_plan);
#use Test::Deep;
use Test::More qw(no_plan);

my $max = 9999;
my $numiter = 10;

##########################################################
# range_before, range_after

for (1 .. $numiter){
    my ($x1,$x2,$y1,$y2) = sort { $a <=> $b } map { int(rand($max)) } (1..4);
    
    my $x = [$x1, $x2]; my $y = [$y1, $y2];
    ok(range_before($x,$y) && range_after($y,$x) && ! range_overlap($x,$y),
        "range_before, range_after"
    );
}

##########################################################
# range_before, range_after

for (1 .. $numiter){
    my ($w,$x,$y,$z) = sort { $a <=> $b } map { int(rand($max)) } (1..4);
    my @overlappers = (
        [[$w,$y],[$x,$z]], # >=1 overlap
        [[$w,$z],[$x,$y]], # containment
        [[$w,$x],[$x,$z]], # single base overlap
        [[$w,$x],[$w,$z]], # left side is same
        [[$w,$y],[$x,$y]], # right side is same
        [[$w,$z],[$w,$z]], # same
    );
    for (@overlappers){
        ok( ! range_before($_->[0],$_->[1]) && 
            ! range_after($_->[1],$_->[0]) && 
            range_overlap($_->[0],$_->[1]) && 
            overlap_ratio($_->[1],$_->[0])>=0 &&
            overlap_ratio($_->[1],$_->[0])<=1 ,
            "range_overlap"
        );
    }

}

##########################################################
# unique range

sub rand_ranges{
    my $len ||= 100;
    return map {
        my $x1 = int(rand($max));
        my $x2 = $x1 + int(rand($max));
        [$x1,$x2];
    } (1 .. $len);
}

sub add_rand_dups{
    my $l = scalar @_;
    for (1 .. int($l/3)){
        my $dup = $_[int(rand($l))];
        push @_, $dup;
    }
    return shuffle @_;
}

sub range_cmp{
    $_[0][0] <=> $_[1][0]
}

for (1 .. $numiter){
    my @original = rand_ranges;
    my @withdupes = add_rand_dups(@original);
    die "add_rand_dups not working..." unless (@original < @withdupes);
    my @uniqed = @{uniq_ranges(\@withdupes)};

    is_deeply([sort \&range_cmp, @original],  [sort \&range_cmp, @uniqed], 
        "make sure uniq_ranges is actually getting rid of duplicates");
}

