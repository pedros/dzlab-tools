#!/usr/bin/perl

use Data::Dumper;
use Getopt::Long;
use strict;
use warnings;
use diagnostics;
disable diagnostics;

# Globals, passed as command line options
my $gfffile = "";
my $readsize = 0;
my $output = "-";
my $verbose = 0;
my $quiet = 0;
my $usage = 0;

# Initial check of command line parameters
&usage;

my $result = GetOptions (
    "gff|f=s" => \$gfffile,
    "readsize|r=i" => \$readsize,
    "output|o:s" => \$output,
    "verbose|v" => sub {enable diagnostics;use warnings;},
    "quiet|q" => sub {disable diagnostics;no warnings;},
    "usage|help|h" => \&usage
    );

# opens gff file
open(my $GFF, "<", $gfffile) or die("Can't read file: $gfffile");

# redirects STDOUT to file if specified by user
if(!($output eq '-')) {
    open(STDOUT, ">", "$output") or die("Can't redirect STDOUT to file: $output");
}

my %HoH=();
while(<$GFF>) {
    chomp($_);
    my %record = %{&readGFF($_)};

    next if($record{'score'} > 1);
    $record{'feature'} =~ m/([ACGTN]{$readsize})/;
    my @methylated = split(//, $1);
    my @unmethylated = split(//, (split "=", $record{'attribute'})[1]);
    
    for(my $i=$record{'start'};$i<$record{'end'};$i++) {
	my $j = $i - $record{'start'};

	if($record{'feature'} =~ m/\/1:[ACGTN]{$readsize}/ && $unmethylated[$j+2] =~ m/[Cc]/) {

	    if($methylated[$j] =~ m/[Cc]/) {
		$HoH{$i}{'c_count'}++;
	    }
	    elsif($methylated[$j] =~ m/[Tt]/) {
		$HoH{$i}{'t_count'}++;
	    }

	    if($unmethylated[$j+3] =~ m/[Gg]/) {
		$HoH{$i}{'cg_count'}++;
	    }
	    elsif(join("",@unmethylated[($j+3)..($j+4)]) =~ m/[^Gg][Gg]/) {
		$HoH{$i}{'chg_count'}++;
	    }
	    elsif(join("",@unmethylated[($j+3)..($j+4)]) =~ m/[^Gg][^Gg]/) {
		$HoH{$i}{'chh_count'}++;
	    }

	    $HoH{$i}{'coord'}=$i;
	    $HoH{$i}{'chr'}=$record{'seqname'};
	    $HoH{$i}{'strand'}=$record{'strand'};
	}
	
# 	elsif($record{'feature'} =~ m/\/2:[ACGTN]{$readsize}/ && $unmethylated[$j+2] =~ m/[Gg]/) {
	    
# 	    if($methylated[$j] =~ m/[Gg]/) {
# 		$HoH{$i}{'c_count'}++;
# 	    }
# 	    elsif($methylated[$j] =~ m/[Aa]/) {
# 		$HoH{$i}{'t_count'}++;
# 	    }

# 	    if($unmethylated[$j+1] =~ m/[Gg]/) {
# 		$HoH{$i}{'cg_count'}++;
# 	    }
# 	    elsif(join("",@unmethylated[$j..($j+1)]) =~ m/[Cc][^Cc]/) {
# 		$HoH{$i}{'chg_count'}++;
# 	    }
# 	    elsif(join("",@unmethylated[$j..($j+1)]) =~ m/[^Cc][^Cc]/) {
# 		$HoH{$i}{'chh_count'}++;
# 	    }

# 	    $HoH{$i}{'coord'}=$i;
# 	    $HoH{$i}{'chr'}=$record{'seqname'};
# 	}
    }
}

for my $i (sort keys %HoH) {
    if(!exists $HoH{$i}{'c_count'} && !exists $HoH{$i}{'t_count'}) {
	$HoH{$i}{'score'} = 0;
    }
    else {
	$HoH{$i}{'score'} = $HoH{$i}{'c_count'}/($HoH{$i}{'c_count'}+$HoH{$i}{'t_count'});
    }

    my $context = ".";
    if($HoH{$i}{'cg_count'}>$HoH{$i}{'chg_count'} && $HoH{$i}{'g_count'}>$HoH{$i}{'chh_count'}) {
	$context = "CG";
    }
    elsif($HoH{$i}{'chh_count'}>$HoH{$i}{'chg_count'}) {
	$context = "CHH";
    }
    elsif($HoH{$i}{'chg_count'}>0) {
	$context = "CHG";
    }

    print join("\t",
	       $HoH{$i}{'chr'},
	       "dz_cm",
	       $context,
	       $HoH{$i}{'coord'},
	       $HoH{$i}{'coord'},
	       sprintf("%.3f", $HoH{$i}{'score'}),
	       $HoH{$i}{'strand'},
	       ".",
	       join(";", "c=$HoH{$i}{'c_count'}", "t=$HoH{$i}{'t_count'}")
	), "\n";
}

#print Dumper(\%HoH);




close($GFF);
close(STDOUT);
# done


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

# sub writeGFF {
#     my %rec = \$_;
#     return
# 	join("\t",
# 	     $rec{'seqname'},
# 	     $rec{'source'},
# 	     $rec{'feature'},
# 	     $rec{'start'},
# 	     $rec{'end'},
# 	     $rec{'score'},
# 	     $rec{'strand'},
# 	     $rec{'frame'},
# 	     $rec{'attribute'}
# 	);
# }

sub usage {
    if ($usage || @ARGV<2) {
	print STDERR
	    "countMethylation.pl <PARAMETERS> [OPTIONS]
\t<--gff>\tGFF alignment input file
\t<--readsize>\tNo. of bps per read
\t[--output]\tFilename to write results to (default is STDOUT)
\t[--verbose]\tOutput perl's diagnostic and warning messages
\t[--quiet]\tSupress perl's diagnostic and warning messages
\t[--usage]\tPrint this information
Takes a GFF-formatted input file with alignment information.
Calculates methylation levels and contexts for each cytosine in scaffold.
Methylated sequence is assumed to be in last field of 'feature' column.
Original reference sequence is assumed to be in the 'attribute' column.
";
	exit 1;
    }
}
