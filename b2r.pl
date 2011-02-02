#!/usr/bin/env perl
# ___UNDOCUMENTED___

use strict;
use warnings;

unless (@ARGV == 1) {
    print STDERR ("Usage: b2r.pl <fastafile.fas>\n") ;
    exit(-1);
}
my $fastafile = "$ARGV[0]";
my @fastaseq=();

open(FAS,"<$fastafile") or die "Can't read input file";
while(my $line = <FAS>) {
    if($line !~ m/^>/) {
	my @chars = split(//, $line);
	if(!int(rand(9)==0)) {
	    	push @fastaseq, $line;   
		next;
	}
	if(int(rand(1)==0)) {
	    $chars[int(rand(length(@chars)))-1]=chr(int(rand(3)) + 65);
	    $line = join("", @chars); 
	    push @fastaseq, $line;
	    next;
	}
	if(int(rand(2)==0)) {
	    $chars[int(rand(length(@chars)))-1]=chr(int(rand(3)) + 65);
	    $chars[int(rand(length(@chars)))-1]=chr(int(rand(3)) + 65);
	    $line = join("", @chars); 
	    push @fastaseq, $line;
	    next;
	}  
	if(int(rand(3)==0)) {
	    $chars[int(rand(length(@chars)))-1]=chr(int(rand(3)) + 65);
	    $chars[int(rand(length(@chars)))-1]=chr(int(rand(3)) + 65);
	    $chars[int(rand(length(@chars)))-1]=chr(int(rand(3)) + 65);
	    $line = join("", @chars); 
	    push @fastaseq, $line;
	    next;
	}	
    }
}

print STDERR chr(int(rand(3)) + 65);

close(FAS);
print @fastaseq;
