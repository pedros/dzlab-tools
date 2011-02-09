#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Test::More qw(no_plan);
use GFF;

do_gff {
    say gff_to_string $_;
} file => 't/test1.gff', locus => 'ID';

ok(1);
