#!/usr/bin/env perl
# ___UNDOCUMENTED___

use strict;
use warnings;

unless (@ARGV == 2) {
    print STDERR ("Usage: convert.pl [c2t|g2a] <fastafile.fas>\n") ;
    exit(-1);
}

my $pattern   = $ARGV[0];
my $fastafile = $ARGV[1];

open(FAS,"<$fastafile") or die "Can't read input file";
while(my $line = <FAS>) {
    if($line =~ m/[ACGTN]+/i) {
	$line =~ tr/Cc/Tt/ if $pattern eq 'c2t';
	$line =~ tr/Cc/Tt/ if $pattern eq 'g2a';
    }
    print $line;
}
close(FAS);
