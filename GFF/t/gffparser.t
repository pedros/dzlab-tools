#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Test::More qw(no_plan);
use GFF::Parser;
use GFF;

my $it = GFF::Parser->new(file => 't/test1.gff',locus => 'ID');

while (defined (my $gff = $it->next())){
    if (is_gff $gff){
        print Dumper $gff;
    }
}

ok(1);
