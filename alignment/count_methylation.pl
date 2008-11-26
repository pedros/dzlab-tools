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
	   "Sequence",
	   "Scaffold",
	   "Start",
	   "End",
	   "Ratio",
	   "Strand",
	   "Context",
	   "C;T;CG;CHG;CHH;TG;THG;THH",
    ), "\n";

while(<$GFF>) {
    chomp($_);
    my %record = %{&readGFF($_)};
    next if($record{'score'} > 1);

    if($record{'score'} > 1) {
    print join("\t",
	       $record{'seqname'},
	       $record{'source'},
	       $record{'feature'},
	       $record{'start'},
	       $record{'end'},
	       $record{'score'},
	       $record{'strand'},
	       $record{'frame'},
	       $record{'attribute'}), "\n";
    next;
    }

    $record{'feature'} =~ m/([ACGTN]{$readsize})/;   
    my @methylated = split(//, $1);
    my @unmethylated = split(//, (split "=", $record{'attribute'})[1]);
    my ($c_count, $t_count, $cg_count, $chg_count, $chh_count, $tg_count, $thg_count, $thh_count, $ratio) = (0,0,0,0,0,0,0,0,0);

    for(my $i=$record{'start'};$i<$record{'end'}+2;$i++) {
	my $j = $i - $record{'start'};
	print STDERR $j, "\n";

	if($j > 1 && $record{'feature'} =~ m/\/1:[ACGTN]{$readsize}/ && $unmethylated[$j] =~ m/[Cc]/) {

	    if($methylated[$j-2] =~ m/[Cc]/) {
		$c_count++;
	    }
	    if($methylated[$j-2] =~ m/[Tt]/) {
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

	    if(join("",@unmethylated[$j..($j+1)]) =~ m/[Tt][Gg]/) {
		$tg_count++;
	    }
	    if(join("",@unmethylated[$j..($j+2)]) =~ m/[Tt][^Gg][Gg]/) {
		$thg_count++;
	    }
	    if(join("",@unmethylated[$j..($j+2)]) =~ m/[Tt][^Gg][^Gg]/) {
		$thh_count++;
	    }

	}

	if($j < $record{'end'} && $record{'feature'} =~ m/\/2:[ACGTN]{$readsize}/ && $unmethylated[$j] =~ m/[Gg]/) {

	    if($j > 1 && $methylated[$j-2] =~ m/[Gg]/) {
		$c_count++;
	    }
	    if($j > 1 && $methylated[$j-2] =~ m/[Aa]/) {
		$t_count++;
	    }

	    if($j > 0 && join("",@unmethylated[$j..($j+1)]) =~ m/[Cc][Gg]/) {
		$cg_count++;
	    }
	    if(join("",@unmethylated[$j..($j+2)]) =~ m/[Cc][^Cc][Gg]/) {
		$chg_count++;
	    }
	    if(join("",@unmethylated[$j..($j+2)]) =~ m/[^Cc][^Cc][Gg]/) {
		$chh_count++;
	    }

	    if($j > 0 && join("",@unmethylated[$j..($j+1)]) =~ m/[Cc][Aa]/) {
		$tg_count++;
	    }
	    if(join("",@unmethylated[$j..($j+2)]) =~ m/[Cc][^Cc][Aa]/) {
		$thg_count++;
	    }
	    if(join("",@unmethylated[$j..($j+2)]) =~ m/[^Cc][^Cc][Aa]/) {
		$thh_count++;
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
    
    if($c_count == 0 && $t_count == 0) {
	$ratio = 0;
    }
    else {
	$ratio = ($c_count/($c_count+$t_count));
    }

    print join("\t",
	       $record{'seqname'},
	       $record{'feature'},
	       join("",@unmethylated),
	       $record{'start'},
	       $record{'end'},
	       sprintf("%.5f", $ratio),
	       $record{'strand'},
	       $context,
	       join(";",
		    "c=$c_count",
		    "t=$t_count",
		    "cg=$cg_count",
		    "chg=$chg_count",
		    "chh=$chh_count",
		    "tg=$tg_count",
		    "thg=$thg_count",
		    "thh=$thh_count"
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
