#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;

use List::Util qw(max);

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

    for(my $i=$record{'start'};$i<$record{'end'};$i++) {
	my $j = $i - $record{'start'};

	if($record{'feature'} =~ m/\/1:[ACGTN]{$readsize}/ && $unmethylated[$j+2] =~ m/[Cc]/) {

	    if($methylated[$j] =~ m/[Cc]/) {
		$c_count++;
		if($unmethylated[$j+3] =~ m/[Gg]/) {
		    $cg_count++;
		}
		elsif(join("",@unmethylated[($j+3)..($j+4)]) =~ m/[^Gg][Gg]/) {
		    $chg_count++;
		}
		elsif(join("",@unmethylated[($j+3)..($j+4)]) =~ m/[^Gg][^Gg]/) {
		    $chh_count++;
		}
	    }
	    elsif($methylated[$j] =~ m/[Tt]/) {
		$t_count++;
		if($unmethylated[$j+3] =~ m/[Gg]/) {
		    $tg_count++;
		}
		elsif(join("",@unmethylated[($j+3)..($j+4)]) =~ m/[^Gg][Gg]/) {
		    $thg_count++;
		}
		elsif(join("",@unmethylated[($j+3)..($j+4)]) =~ m/[^Gg][^Gg]/) {
		    $thh_count++;
		}
	    }
	}

	elsif($record{'feature'} =~ m/\/2:[ACGTN]{$readsize}/ && $unmethylated[$j+2] =~ m/[Gg]/) {

	    if($methylated[$j] =~ m/[Gg]/) {
		$c_count++;
		if($unmethylated[$j+1] =~ m/[Cc]/) {
		    $cg_count++;
		}
		elsif(join("",@unmethylated[$j..($j+1)]) =~ m/[Cc][^Cc]/) {
		    $chg_count++;
		}
		elsif(join("",@unmethylated[$j..($j+1)]) =~ m/[^Cc][^Cc]/) {
		    $chh_count++;
		}
		
	    }
	    elsif($methylated[$j] =~ m/[Aa]/) {
		$t_count++;
		if($unmethylated[$j+1] =~ m/[Cc]/) {
		    $tg_count++;
		}
		elsif(join("",@unmethylated[$j..($j+1)]) =~ m/[Cc][^Cc]/) {
		    $thg_count++;
		}
		elsif(join("",@unmethylated[$j..($j+1)]) =~ m/[^Cc][^Cc]/) {
		    $thh_count++;
		}
	    }    
	}
    }
    
    my %contexts=();
    $contexts{'cg_count'}=$cg_count;
    $contexts{'chg_count'}=$chg_count;
    $contexts{'chh_count'}=$chh_count;
    $contexts{'tg_count'}=$tg_count;
    $contexts{'thg_count'}=$thg_count;
    $contexts{'thh_count'}=$thh_count;

    my $maxval = max(values %contexts);
    %contexts = reverse %contexts;
    
    my $context = $contexts{$maxval};
    $context =~ s/_count//;
    $context =~ tr/a-z/A-Z/;


#     my ($context, $c_context, $t_context);
#     if($cg_count>$chg_count && $cg_count>$chh_count) {$c_context = "cg";}
#     elsif($chh_count>$chg_count) {$c_context = "chh";}
#     elsif($chg_count>0) {$c_context = "chg";}

#     if($cg_count>$chg_count && $cg_count>$chh_count) {$t_context = "tg";}
#     elsif($chh_count>$chg_count) {$t_context = "thh";}
#     elsif($chg_count>0) {$t_context = "thg";}

#     if("${c_context}_count" >  "${t_context}_count") {$context = $c_context;}
#     else {$context = $t_count;}
    
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
	       $context ,
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
