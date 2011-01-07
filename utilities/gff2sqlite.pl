#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';

use FindBin;
use lib "$FindBin::Bin/../DZLab-Tools/lib";
use DZLab::Tools::GFF qw/gff_make_iterator/;
use DZLab::Tools::GFFStore;

my $gffstore = DZLab::Tools::GFFStore->new({
        attributes => {ID => 'text', Parent => 'text'},
        verbose => 1,
    });

$gffstore->slurp({handle => \*ARGV});

