#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;

use List::Util qw(max);

die("Usage:\tcount_methylation.pl <gff alignment file> <read size>") unless @ARGV==2;

open(my $GFF, "<", $ARGV[0]);
my $readsize = $ARGV[1];

print join("\t",
	   "#Chr",
	   "Original",
	   "Methylated",
	   "Start",
	   "End",
	   "Ratio",
	   "Strand",
	   "Context",
	   "C;T;CG;CHG;CHH",
    ), "\n";

while(<$GFF>) {
    chomp($_);
    my %record = %{&readGFF($_)};
    next unless($record{'score'} == 1);
    $record{'feature'} =~ m/([ACGT]{$readsize})/;   
    my @methylated = split(//, $1);
    my @unmethylated = split(//, (split "=", $record{'attribute'})[1]);
    my ($c_count, $t_count, $cg_count, $chg_count, $chh_count) = (0,0,0,0,0);

    for(my $i=$record{'start'};$i<$record{'end'};$i++) {
	my $j = $i - $record{'start'};
	if($unmethylated[$j] =~ m/[CcGg]/) {
	    if($record{'strand'} == '+') {
		if($methylated[$j] =~ m/[Cc]/) {
		    $c_count++;
		}
		if($methylated[$j] =~ m/[Tt]/) {
		    $t_count++;
		}
		if(join("",@unmethylated[$j..($j+1)]) =~ m/[Cc][Gg]/) {
		    $cg_count++;
		}
		if(join("",@unmethylated[$j..($j+2)]) =~ m/[Cc][^Gg][Gg]/) {
		    $chg_count++;
		}
		if(join("",@unmethylated[$j..($j+2)]) =~ m/[Cc][^Gg][^Gg]/) {
		    $chh_count++;
		}
	    }

	    if($record{'strand'} == '-') {
		if($methylated[$j] =~ m/[Gg]/) {
		    $c_count++;
		}
		if($methylated[$j] =~ m/[Aa]/) {
		    $t_count++;
		}
		if(join("",@unmethylated[$j..($j+1)]) =~ m/[Cc][Gg]/) {
		    $cg_count++;
		}
		if(join("",@unmethylated[$j..($j+2)]) =~ m/[Cc][^Gg][Gg]/) {
		    $chg_count++;
		}
		if(join("",@unmethylated[$j..($j+2)]) =~ m/[^Cc][^Cc][Gg]/) {
		    $chh_count++;
		}
	    }
	}	   
    }
    
    my $context = ".";
    if($cg_count>$chg_count && $cg_count>$chh_count) {
	$context = "CG";
    }
    elsif($chh_count>$chg_count) {
	$context = "CHH";
    }
    elsif($chg_count>0) {
	$context = "CHG";
    }
    
    print join("\t",
	       $record{'seqname'},
	       $record{'feature'},
	       join("",@unmethylated),
	       $record{'start'},
	       $record{'end'},
	       ($c_count/($c_count+$t_count)),
	       $record{'strand'},
	       $context,
	       join(";",
		    "c=$c_count",
		    "t=$t_count",
		    "cg=$cg_count",
		    "chg=$chg_count",
		    "chh=$chh_count"
	       ),
	       ), "\n";
}


sub readGFF {
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $_[0]);
    my %rec = (
	'seqname'=>$seqname,
	'source'=>$source,
	'feature'=>$feature,
	'start'=>$start,
	'end'=>$end,
	'score'=>$score,
	'strand'=>$strand,
	'frame'=>$strand,
	'attribute'=>$attribute
	);
    return \%rec;
}
