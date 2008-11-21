#!/usr/bin/env perl

use Data::Dumper;
use Getopt::Long;
use strict;
use warnings;
use diagnostics;
#disable diagnostics;
#no warnings;

# Globals, passed as command line options
my $gfffile = "";
my $width = 0;
my $step = 0;
my $readsize = 0;
my $output = "-";
my $verbose = 0;
my $quiet = 0;
my $usage = 0;

# Initial check of command line parameters
&usage;

my $result = GetOptions (
    "gff|f=s" => \$gfffile,
    "width|w=i" => \$width,
    "step|s=i" => \$step,
    "readsize|r=i" => \$readsize,
    "output|o:s" => \$output,
    "verbose|v" => sub {enable diagnostics;use warnings;},
    "quiet|q" => sub {no warnings;disable diagnostics;},
    "usage|help|h" => \&usage
    );

# opens gff file
open(my $GFF, "<", $gfffile) or die("Can't read file: $gfffile");

# redirects STDOUT to file if specified by user
if(!($output eq '-')) {
    open(STDOUT, ">", "$output") or die("Can't redirect STDOUT to file: $output");
}

while(<$GFF>) {
    chomp($_);
    my %record = %{&readGFF($_)};

    if($record{'score'}!=1){next;}

    $record{'feature'} =~ m/([ACGT]){$readsize}/;
    
    my @methylated = split(//, $1);
    my @unmethylated = split(//, (split "=", $record{'attribute'})[1]);
    
    my %HoH=();

    for(my $i=$record{'start'};$i<$record{'end'};$i++) {
	if($unmethylated[$i] =~ m/[CcGg]/) {
	    if($record{'strand'} == '+') {
		if($methylated[$i] =~ m/[Cc]/) {
		    $HoH{$i}{'c_count'}++;
		}
		if($methylated[$i] =~ m/[Tt]/) {
		    $HoH{$i}{'t_count'}++;
		}
		if(join(@methylated[$i..($i+1)]) =~ m/[Cc][Gg]/) {
		    $HoH{$i}{'cg_count'}++;
		}
		if(join(@methylated[$i..($i+2)]) =~ m/[Cc][^Gg][Gg]/) {
		    $HoH{$i}{'chg_count'}++;
		}
		if(join(@methylated[$i..($i+2)]) =~ m/[Cc][^Gg][^Gg]/) {
		    $HoH{$i}{'chh_count'}++;
		}
		print $HoH{$i}{'c_count'}++;
	    }

	    if($record{'strand'} == '-') {
		if($methylated[$i] =~ m/[Gg]/) {
		    $HoH{$i}{'c_count'}++;
		}
		if($methylated[$i] =~ m/[Aa]/) {
		    $HoH{$i}{'t_count'}++;
		}
		if(join(@methylated[$i..($i+1)]) =~ m/[Cc][Gg]/) {
		    $HoH{$i}{'cg_count'}++;
		}
		if(join(@methylated[$i..($i+2)]) =~ m/[Cc][^Gg][Gg]/) {
		    $HoH{$i}{'chg_count'}++;
		}
		if(join(@methylated[$i..($i+2)]) =~ m/[^Cc][^Cc][Gg]/) {
		    $HoH{$i}{'chh_count'}++;
		}
	    }
	}	   
    }  
}

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
    if ($usage || @ARGV<4) {
	print STDERR
	    "countMethylation.pl <PARAMETERS> [OPTIONS]
\t<--gff>\tGFF alignment input file
\t<--width>\tSliding window width
\t<--step>\tSlidting window step size
\t<--readsize>\tNo. of bps per read
\t[--output]\tFilename to write results to (default is STDOUT)
\t[--verbose]\tOutput perl's diagnostic and warning messages
\t[--quiet]\tSupress perl's diagnostic and warning messages
\t[--usage]\tPrint this information
Takes a GFF-formatted input file with alignment information and calculates methylation levels.
Methylated sequence is assumed to be in last field of 'feature' column.
Original reference sequence is assumed to be in the 'attribute' column.
";
	exit 1;
    }
}
