#!/usr/bin/env perl

use strict;
use warnings;
use Test::Simple qw(no_plan);
#use Test::More qw(no_plan);
use List::Util 'shuffle';

my $orig = 't/sorted.gff';
my $dirty = 't/unsorted.gff';
my $test = 't/test.gff';


for (1..10){
    shuffle_file($orig,$dirty);
    system("./sort_gff.pl -i $dirty -o $test");
    ok(files_identical($test,$orig));
}

unlink $test;
unlink $dirty;

sub shuffle_file{
    my ($in,$out) = @_;
    open my $inh, '<', $in;
    open my $outh, '>', $out;
    print $outh shuffle(<$inh>);
    close $outh;
    close $inh;
}

sub files_identical{
    local $/=undef;
    my ($f1,$f2) = @_;
    open my $f1h, '<', $f1;
    open my $f2h, '<', $f2;
    my $s1 = <$f1h>;
    my $s2 = <$f2h>;
    close $f1h;
    close $f2h;
    return $s1 eq $s2;
}
