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

# Grabs and parses command line options
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

# prints out a commented line with header fields
print join("\t",
	   "#Chr",
	   "Source",
	   "Context",
	   "Start",
	   "End",
	   "Ratio",
	   "Strand",
	   "Frame",
	   "C;T;CG;CHG;CHH;TG;THG;THH",
    ), "\n";

# declare and initialize a hash of hashes, where each inner hash's key is a 'C' coordinate
# and its key/value pairs are total c/t counts, and contexts
my %HoH=();

# reads in input file until EOF
while(<$GFF>) {
    chomp($_);

    # parses each record line into a simple hash
    my %record = %{&readGFF($_)};

    # skips records with non-unique matches
    next if($record{'score'} != 1);

    # grabs methylated sequence from gff file
    # assumes sequence is exactly $readsize bps long
    # note that we can do this because the scaffold read is 4 bps longer
    $record{'feature'} =~ m/([ACGTN]{$readsize})/;
    my @methylated = split(//, $1);
    # grabs scaffold sequence from gff file
    # assumes sequence is in the 'attribute' field, separated from extraneous information by a '='
    my @unmethylated = split(//, (split "=", $record{'attribute'})[1]);
    
    # loops through each character in current sequences
    for(my $i=$record{'start'};$i<$record{'end'};$i++) {
	# sets $j to 0 for forward strand coordinates
	my $j = $i - $record{'start'};
	# sets $k to $record{'end'} for reverse strand coordinates
	my $k = $record{'end'} - $j;

	my $coord = $i;
	$coord = $k if($record{'strand'} eq '-');

	# checks that we're looking at a left/1 sequence AND current character is a 'C'
	# this regex is adapted to the Solexa sequences we have, so it might change in the future
	# we're looking 2 characters ahead in the scaffold because it has 2 extra chars in the beginning
 	if($unmethylated[$j+2] =~ m/[Cc]/ && 
 	   $record{'feature'} =~ m/\/1:[ACGTN]{$readsize}/) {

	    # checks what happened in the methylated sequence
	    # and updates the appropriate 'c' or 't' count
	    if($methylated[$j] =~ m/[Cc]/) {
		$HoH{$coord}{'c_count'}++;
	    }
	    elsif($methylated[$j] =~ m/[Tt]/) {
		$HoH{$coord}{'t_count'}++;
	    }

	    # checks the context by looking ahead +3 or +3..+4 bps
	    # because the scaffold read is displaced by 2 bps
	    if($unmethylated[$j+3] =~ m/[Gg]/) {
		$HoH{$coord}{'cg_count'}++;
	    }
	    elsif(join("",@unmethylated[($j+3)..($j+4)]) =~ m/[^Gg][Gg]/) {
		$HoH{$coord}{'chg_count'}++;
	    }
	    elsif(join("",@unmethylated[($j+3)..($j+4)]) =~ m/[^Gg][^Gg]/) {
		$HoH{$coord}{'chh_count'}++;
	    }

	    # grab some necessary information into the data structure
	    # we will print this information later
	    $HoH{$coord}{'coord'} = $coord;
	    $HoH{$coord}{'chr'} = $record{'seqname'};
	    $HoH{$coord}{'strand'} .= $record{'strand'};
	    $HoH{$coord}{'end'} .= "/1";	
	}
	
	# checks that we're looking at a right/2 sequence AND current character is a 'C'
	# this regex is adapted to the Solexa sequences we have, so it might change in the future
	# we're looking 2 characters before in the scaffold because it has 2 extra chars in the beginning
	# AND we're looking in the reverse strand, so the context goes right-to-left
	elsif($unmethylated[$j+2] =~ m/[Gg]/ && 
	      $record{'feature'} =~ m/\/2:[ACGTN]{$readsize}/) {

	    # checks what happened in the methylated sequence
	    # and updates the appropriate 'c' or 't' count
	    if($methylated[$j] =~ m/[Gg]/) {
		$HoH{$coord}{'c_count'}++;
	    }
	    elsif($methylated[$j] =~ m/[Aa]/) {
		$HoH{$coord}{'t_count'}++;
	    }

	    # checks the context by looking behind +1 or +0..+1 bps
	    # because the scaffold read is displaced by 2 bps
	    if($unmethylated[$j+1] =~ m/[Cc]/) {
		$HoH{$coord}{'cg_count'}++;
	    }
	    elsif(join("",@unmethylated[$j..($j+1)]) =~ m/[Cc][^Cc]/) {
		$HoH{$coord}{'chg_count'}++;
	    }
	    elsif(join("",@unmethylated[$j..($j+1)]) =~ m/[^Cc][^Cc]/) {
		$HoH{$coord}{'chh_count'}++;
	    }

	    # grab some necessary information into the data structure
	    # we will print this information later
	    $HoH{$coord}{'coord'} = $coord;
	    $HoH{$coord}{'chr'} = $record{'seqname'};
	    $HoH{$coord}{'strand'} .= $record{'strand'};
	    $HoH{$coord}{'end'} .= "/2";
	}
    }
}

# loops through every initialized key in main hash
# keys are sorted, so the output is going to be sorted by starting coordinate
for my $i (sort {$a <=> $b} keys %HoH) {

    # we need to check that a given key/value pair was initialized
    # if it wasn't, initialize it to zero (to avoid division-by-zero, etc)
    if(!exists $HoH{$i}{'c_count'}) {$HoH{$i}{'c_count'} = 0;}
    if(!exists $HoH{$i}{'t_count'}) {$HoH{$i}{'t_count'} = 0;}
    if(!exists $HoH{$i}{'cg_count'}) {$HoH{$i}{'cg_count'} = 0;}
    if(!exists $HoH{$i}{'chg_count'}) {$HoH{$i}{'chg_count'} = 0;}
    if(!exists $HoH{$i}{'chh_count'}) {$HoH{$i}{'chh_count'} = 0;}

    # calculates c/(c+t), which varies between 0 and 1
    # 0 is no methylation, 1 is total methylation
    if($HoH{$i}{'c_count'} + $HoH{$i}{'t_count'} == 0) {$HoH{$i}{'score'} = 0;}
    else {$HoH{$i}{'score'} = $HoH{$i}{'c_count'}/($HoH{$i}{'c_count'}+$HoH{$i}{'t_count'});}

    # finds the appropriate context to put in 'feature' field
    my $context = ".";
    if($HoH{$i}{'cg_count'}>$HoH{$i}{'chg_count'} && $HoH{$i}{'cg_count'}>$HoH{$i}{'chh_count'}) {$context = "CG";}
    elsif($HoH{$i}{'chh_count'}>$HoH{$i}{'chg_count'}) {$context = "CHH";}
    elsif($HoH{$i}{'chg_count'}>0) {$context = "CHG";}

    if(($HoH{$i}{'cg_count'} > 0 && $HoH{$i}{'chg_count'} > 0) ||
       ($HoH{$i}{'cg_count'} > 0 && $HoH{$i}{'chh_count'} > 0) ||
       ($HoH{$i}{'chg_count'} > 0 && $HoH{$i}{'chh_count'} > 0)) {$context=$HoH{$i}{'end'};}

    # prints a single gff record
    print join("\t",
	       $HoH{$i}{'chr'},
	       "dz_cm",
	       $context,
	       $HoH{$i}{'coord'},
	       $HoH{$i}{'coord'},
	       sprintf("%.3f", $HoH{$i}{'score'}),
	       $HoH{$i}{'strand'},
	       ".",
	       join(";", "c=$HoH{$i}{'c_count'}", "t=$HoH{$i}{'t_count'}", "cg=$HoH{$i}{'cg_count'}", "chg=$HoH{$i}{'chg_count'}", "chh=$HoH{$i}{'chh_count'}")
	), "\n";

    delete $HoH{$i};
}

# debugging information
#print Dumper(\%HoH);

close($GFF);
close(STDOUT);
# done

# readGFF reads in a single GFF line, and outputs a hash reference
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


# prints out usage information
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
