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
use GFF::Sort;
use GFF::Parser;

my $orig = 't/sorted.gff';
my $dirty = 't/unsorted.gff';
my $test = 't/unsorted.gff.sorted';


for (1..10){
    shuffle_file($orig,$dirty);
    gff_sort(file => $dirty, overwrite => 0, cols => ['start'], manual => 0);
    #system("./sort_gff.pl -i $dirty -o $test");
    files_identical($test,$orig);
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
    my ($f1,$f2) = @_;
    my ($gff1, $gff2) = map { GFF::Parser->new(file => $_)->slurp() } ($f1, $f2);
    is_deeply($gff1,$gff2,"$f1 and $f2 identical contents");
}
