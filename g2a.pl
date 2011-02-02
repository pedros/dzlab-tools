#!/usr/bin/env perl
# ___UNDOCUMENTED___

use strict;
use warnings;

unless (@ARGV == 1) {
    print STDERR ("Usage: g2a.pl <fastafile.fas>\n") ;
    exit(-1);
}
my $fastafile = "$ARGV[0]";

open(FAS,"<$fastafile") or die "Can't read input file";
while(my $line = <FAS>) {
    if($line !~ m/^>/) {
	$line =~ tr/Gg/Aa/;
    }
    print $line;
}
close(FAS);
