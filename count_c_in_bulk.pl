#!/usr/bin/env perl
# ___UNDOCUMENTED___

use strict;

if (@ARGV<3) {
    die("Usage: countC.pl <gff alignment file> <type (exclude, include>) <regex to include or exclude>\n");
}

my $gfffile=$ARGV[0];
my $type=$ARGV[1];
my $regex=$ARGV[2];

open(my $GFF, "<", "$gfffile") or die("Can't open $gfffile");

my ($total_read1_C, $total_read1_CG, $total_read1_CHG, $total_read1_CHH, $total_read1_BP, $independent_read1_C, $independent_read1_T, $independent_read1_A, $independent_read1_G, $independent_total_read1) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
my ($total_read2_C, $total_read2_CG, $total_read2_CHG, $total_read2_CHH, $total_read2_BP, $independent_read2_C, $independent_read2_T, $independent_read2_A, $independent_read2_G, $independent_total_read2) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
my ($total_gene_C, $total_gene_CG, $total_gene_CHG, $total_gene_CHH, $total_gene_BP, $independent_gene_C, $independent_gene_T, $independent_gene_A, $independent_gene_G, $independent_total_gene) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

while (my $gffline = <$GFF>) {

    next if ($gffline =~ m/^#.*$|^\s*$/);
    next if ((split "\t", $gffline)[5]!=1);
    next if ((split "\t", $gffline)[0]!~m/$regex/i && $type=~m/include/i);
    next if ((split "\t", $gffline)[0]=~m/$regex/i && $type=~m/exclude/i);

    chomp $gffline;

    my $readline = (split "\t", $gffline)[2];
    my $geneline = (split "\t", $gffline)[8];

    my $readnumber;
    ($readnumber, $readline) = $readline =~ m/\/([0-9]):([ACGTN]+)/; 
    $geneline = (split '=', $geneline)[1]; 

    my @readseq=split(//, $readline);
    my @geneseq=split(//, $geneline);

    $total_read1_BP = $total_read1_BP + scalar(@readseq) if ($readnumber == 1);
    $total_read2_BP = $total_read2_BP + scalar(@readseq) if ($readnumber == 2);
    $total_gene_BP = $total_gene_BP + scalar(@geneseq) - 4;
    
    for (my $i=0;$i<@readseq;$i++) {

        if ($readnumber == 1) {
            $independent_total_read1++;
        } else {
            $independent_total_read2++;
        }
        
	if ( ($readseq[$i] =~ m/[Cc]/ && $readnumber == 1) ||
                 ($readseq[$i] =~ m/[Gg]/ && $readnumber == 2) ) {
	    $independent_read1_C++ if ($readnumber == 1);
	    $independent_read2_C++ if ($readnumber == 2);
	} elsif ( ($readseq[$i] =~ m/[Tt]/ && $readnumber == 1) ||
                      ($readseq[$i] =~ m/[Aa]/ && $readnumber == 2) ) {
	    $independent_read1_T++ if ($readnumber == 1);
	    $independent_read2_T++ if ($readnumber == 2);
	} elsif ( ($readseq[$i] =~ m/[Gg]/ && $readnumber == 1) ||
                      ($readseq[$i] =~ m/[Cc]/ && $readnumber == 2) ) {
	    $independent_read1_G++ if ($readnumber == 1);
	    $independent_read2_G++ if ($readnumber == 2);
	} elsif ( ($readseq[$i] =~ m/[Aa]/ && $readnumber == 1) ||
                      ($readseq[$i] =~ m/[Tt]/ && $readnumber == 2) ) {
	    $independent_read1_A++ if ($readnumber == 1);
	    $independent_read2_A++ if ($readnumber == 2);
	}
        
	if ((join("",@readseq[$i..($i+1)])=~m/[Cc][Gg]/ && $readnumber == 1) ||
                (join("",@readseq[($i-1)..$i])=~m/[Cc][Gg]/ && $readnumber == 2)) {
	    if ($readnumber == 1) {
		$total_read1_C++;
		$total_read1_CG++;
	    } else {
		$total_read2_C++;
		$total_read2_CG++;
	    }		
	} elsif ((join("",@readseq[$i..($i+2)])=~m/[Cc][^Gg][Gg]/ && $readnumber == 1) ||
                     (join("",@readseq[($i-2)..$i])=~m/[Cc][^Cc][Gg]/ && $readnumber == 2)) {
	    if ($readnumber == 1) {
		$total_read1_C++;
		$total_read1_CHG++;
	    } else {
		$total_read2_C++;
		$total_read2_CHG++;	
            }		
	} elsif ((join("",@readseq[$i..($i+2)])=~m/[Cc][^Gg][^Gg]/ && $readnumber == 1) ||
                     (join("",@readseq[($i-2)..$i])=~m/[^Cc][^Cc][Gg]/ && $readnumber == 2)) {
	    if ($readnumber == 1) {
		$total_read1_C++;
		$total_read1_CHH++;
	    } else {
		$total_read2_C++;
		$total_read2_CHH++;
	    }		
	}
    }

    for (my $i=2;$i<@geneseq-2;$i++) {

        $independent_total_gene++;
        
	if ( ($geneseq[$i] =~ m/[Cc]/ && $readnumber == 1) ||
                 ($geneseq[$i] =~ m/[Gg]/ && $readnumber == 2) ) {
	    $independent_gene_C++ if ($readnumber == 1);
	    $independent_gene_C++ if ($readnumber == 2);
	} elsif ( ($geneseq[$i] =~ m/[Aa]/ && $readnumber == 1) ||
                      ($geneseq[$i] =~ m/[Tt]/ && $readnumber == 2) ) {
	    $independent_gene_A++ if ($readnumber == 1);
	    $independent_gene_A++ if ($readnumber == 2);
	} elsif ( ($geneseq[$i] =~ m/[Gg]/ && $readnumber == 1) ||
                      ($geneseq[$i] =~ m/[Cc]/ && $readnumber == 2) ) {
	    $independent_gene_G++ if ($readnumber == 1);
	    $independent_gene_G++ if ($readnumber == 2);
	} elsif ( ($geneseq[$i] =~ m/[Tt]/ && $readnumber == 1) ||
                      ($geneseq[$i] =~ m/[Aa]/ && $readnumber == 2) ) {
	    $independent_gene_T++ if ($readnumber == 1);
	    $independent_gene_T++ if ($readnumber == 2);
	}
        
	if ((join("",@geneseq[$i..($i+1)])=~m/[Cc][Gg]/ && $readnumber == 1) ||
                (join("",@geneseq[($i-1)..$i])=~m/[Cc][Gg]/ && $readnumber == 2)) {
	    $total_gene_C++;
	    $total_gene_CG++;
	} elsif ((join("",@geneseq[$i..($i+2)])=~m/[Cc][^Gg][Gg]/ && $readnumber == 1) ||
                     (join("",@geneseq[($i-2)..$i])=~m/[Cc][^Cc][Gg]/ && $readnumber == 2)) {
	    $total_gene_C++;
	    $total_gene_CHG++;
	} elsif ((join("",@geneseq[$i..($i+2)])=~m/[Cc][^Gg][^Gg]/ && $readnumber == 1) ||
                     (join("",@geneseq[($i-2)..$i])=~m/[^Cc][^Cc][Gg]/ && $readnumber == 2)) {
	    $total_gene_C++;
	    $total_gene_CHH++;
	}
    }
}

close($GFF);

print "================================================================================\nPost correlation frequencies";
print "\nProcessed file: $ARGV[0]";
print "\nType of count: $type $regex";

print "\n\nTotal read1 BP: $total_read1_BP";
print "\nIndependent total read1: $independent_total_read1";
print "\nIndependent read1 A: $independent_read1_A";
print "\nIndependent read1 G: $independent_read1_G";
print "\nIndependent read1 T: $independent_read1_T";
print "\nIndependent read1 C: $independent_read1_C";
print "\nIndependent read1 ACGT: ", $independent_read1_A + $independent_read1_C + $independent_read1_G + $independent_read1_T;
print "\nTotal read1 C: $total_read1_C";
print "\nTotal read1 CG: $total_read1_CG";
print "\nTotal read1 CHG: $total_read1_CHG";
print "\nTotal read1 CHH: $total_read1_CHH";

print "\n\nTotal read2 BP: $total_read2_BP";
print "\nIndependent total read2: $independent_total_read2";
print "\nIndependent read2 A: $independent_read2_A";
print "\nIndependent read2 G: $independent_read2_G";
print "\nIndependent read2 T: $independent_read2_T";
print "\nIndependent read2 C: $independent_read2_C";
print "\nIndependent read2 C: $independent_read2_C";
print "\nIndependent read2 ACGT: ", $independent_read2_A + $independent_read2_C + $independent_read2_G + $independent_read2_T;
print "\nTotal read2 C: $total_read2_C";
print "\nTotal read2 CG: $total_read2_CG";
print "\nTotal read2 CHG: $total_read2_CHG";
print "\nTotal read2 CHH: $total_read2_CHH";

print "\n\nTotal gene BP: $total_gene_BP";
print "\nIndependent total gene: $independent_total_gene";
print "\nIndependent gene A: $independent_gene_A";
print "\nIndependent gene G: $independent_gene_G";
print "\nIndependent gene T: $independent_gene_T";
print "\nIndependent gene C: $independent_gene_C";
print "\nIndependent gene ACGT: ", $independent_gene_A + $independent_gene_C + $independent_gene_G + $independent_gene_T;
print "\nTotal gene C: $total_gene_C";
print "\nTotal gene CG: $total_gene_CG";
print "\nTotal gene CHG: $total_gene_CHG";
print "\nTotal gene CHH: $total_gene_CHH";

print "\n\nTotal C ratio: ", ($total_read1_C + $total_read2_C) / $total_gene_C;
print "\nTotal CG ratio: ", ($total_read1_CG + $total_read2_CG) / $total_gene_CG;
print "\nTotal CHG ratio: ", ($total_read1_CHG + $total_read2_CHG) / $total_gene_CHG;
print "\nTotal CHH ratio: ", ($total_read1_CHH + $total_read2_CHH) / $total_gene_CHH, "\n================================================================================\n";
