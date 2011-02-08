#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Test::More qw(no_plan);
use GFF::Parser;
use GFF;

my $it = GFF::Parser->new(file => 't/test1.gff',locus => 'ID');
my $seq_index = $it->slurp_index('seqname');

#while (my $gff = $it->next()){
#print Dumper $gff;
#}

#print Dumper $seq_index;

ok(1);
