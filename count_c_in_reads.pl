#!/usr/bin/env perl
# ___UNDOCUMENTED___

use strict;

if(@ARGV<5) {die("Usage: countC.pl <elandfile> <fastafile> <strand (0=forward, 1=reverse)> <type (exclude, include)> <regex to exclude or include>\n");}

my $elandfile=$ARGV[0];
my $fastafile=$ARGV[1];
my $reverse=$ARGV[2];
my $type=$ARGV[3];
my $regex=$ARGV[4];

print STDERR "Processing on the reverse strand...\n" if $reverse;

open(my $ELANDFILE, "<", "$elandfile") or die("Can't open $elandfile");
open(my $FASTAFILE, "<", "$fastafile") or die("Can't open $fastafile");

my ($total_C, $total_CG, $total_CHG, $total_CHH, $total_BP)=(0,0,0,0,0);

my $independent_C = 0;

while(my $elandline = <$ELANDFILE>) {

    my $fastaline = <$FASTAFILE>;
    while($fastaline =~ m/^>|^#.*$|^\s*$/) {
	$fastaline = <$FASTAFILE>;
    }

    next if((split "\t", $elandline)[2] =~ m/NM/);
    next if((split "\t", $elandline)[3]!~m/$regex/i && $type=~m/include/i);
    next if((split "\t", $elandline)[3]=~m/$regex/i && $type=~m/exclude/i);

    chomp $fastaline;

    my @seq=split(//, $fastaline);
    $total_BP = $total_BP + scalar(@seq);
    for(my $i=0;$i<@seq;$i++) {

	$independent_C++ if( ($seq[$i] =~ m/[Cc]/ && !$reverse) ||
			     ($seq[$i] =~ m/[Gg]/ && $reverse) );

	if((join("",@seq[$i..($i+1)])=~m/[Cc][Gg]/ && !$reverse) ||
	   (join("",@seq[($i-1)..$i])=~m/[Cc][Gg]/ && $reverse)) {
	    $total_C++;
	    $total_CG++;
	}
	elsif((join("",@seq[$i..($i+2)])=~m/[Cc][^Gg][Gg]/ && !$reverse) ||
	      (join("",@seq[($i-2)..$i])=~m/[Cc][^Cc][Gg]/ && $reverse)) {
	    $total_C++;
	    $total_CHG++;
	}	
	elsif((join("",@seq[$i..($i+2)])=~m/[Cc][^Gg][^Gg]/ && !$reverse) ||
	      (join("",@seq[($i-2)..$i])=~m/[^Cc][^Cc][Gg]/ && $reverse)) {
	    $total_C++;
	    $total_CHH++;
	}
    }
}
close($ELANDFILE);
close($FASTAFILE);

print "Post alignment frequencies";
print "\nProcessed files: $ARGV[0] and $ARGV[1]";
print "\nType of count: $type $regex";
print "\nTotal BP: $total_BP";
print "\nIndependent C: $independent_C";
print "\nTotal C: $total_C";
print "\nTotal CG: $total_CG";
print "\nTotal CHG: $total_CHG";
print "\nTotal CHH: $total_CHH \n\n";
