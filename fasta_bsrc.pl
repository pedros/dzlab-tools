#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use Fasta qw/bisulfite_convert slurp_fasta/;

my $file = shift // die "usage $0 fasta.txt > fasta.txt.bsrc";


bisulfite_convert($file);


