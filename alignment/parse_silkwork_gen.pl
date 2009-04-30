#!/usr/bin/perl

use Data::Dumper;
use strict;
use warnings;
use diagnostics;

open my $MANIFEST, "<", $ARGV[0];
open my $GENOME, "<", $ARGV[1];


my %manifest = ();
while(<$MANIFEST>) {
    chomp $_;
    next if($_ =~ m/^\s*#/ or $_ =~ m/^\s*$/);
    my ($chr, $scaff, undef, undef) = split(/\t/, $_);
    $manifest{$scaff} = $chr;
}

my %sequences = ();
my $chr;
while(my $line = <$GENOME>) {
    if($line =~ m/^>([^\s]+)\s/) {
	if(exists $manifest{$1}) {
	    $chr = $manifest{$1};
	}
	else {
	    $chr = 0;
	}
	next;
    }
    if($chr) {
	$sequences{$chr} .= $line;
    }
}


foreach my $i (sort keys %sequences) {
    next if($i =~ m/^\s*$/);
    print ">", $i, "\n";
    print $sequences{$i};
}
