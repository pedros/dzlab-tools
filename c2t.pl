#!/usr/bin/perl

use strict;
use warnings;

unless (@ARGV == 1) {
    print STDERR ("Usage: c2t.pl <fastafile.fas>\n") ;
    exit(-1);
}
my $fastafile = "$ARGV[0]";
my @fastaseq=();

open(FAS,"<$fastafile") or die "Can't read input file";
while(my $line = <FAS>) {
    if($line !~ m/^>/) {
	$line =~ tr/Cc/Tt/;
    }
    push @fastaseq, $line;
}
close(FAS);
print @fastaseq;
