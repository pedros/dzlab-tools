#!/usr/bin/perl

use warnings;
use strict;
use diagnostics;

my $loriginal=$ARGV[0];
my $roriginal=$ARGV[1];
my $mutated=$ARGV[2];
my $readsize=$ARGV[3];

open(my $LORIG, "<", "$loriginal") or die("Can't open $loriginal");
open(my $RORIG, "<", "$roriginal") or die("Can't open $roriginal");
open(my $MUTATED, "<", "$mutated") or die("Can't open $mutated");

while(my $mutline = <$MUTATED>) {
    my $lorigline=<$LORIG>;
    while($lorigline=~m/^>/) {
	$lorigline=<$LORIG>;
    }
    my $rorigline=<$RORIG>;
    while($rorigline=~m/^>/) {
	$rorigline=<$RORIG>;
    }

    chomp($lorigline);
    chomp($rorigline);

    $mutline=~s/[ACGT]{$readsize}/$lorigline/ if($mutline=~m/1:\t[ACGT]{$readsize}/);
    $mutline=~s/[ACGT]{$readsize}/$rorigline/ if($mutline=~m/2:\t[ACGT]{$readsize}/);
    
    print $mutline;
}

close($LORIG);
close($RORIG);
close($MUTATED);
