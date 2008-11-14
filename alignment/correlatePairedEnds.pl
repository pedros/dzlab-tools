#!/usr/bin/perl
#
# correlatePairedEnds.pl < --left filename [leftendfile.eland3] >           # 5' raw sequence file
#                                             < --right filename [rightendfile.eland3] >      # 3' raw sequence file
#                                             < --reference [referencefile.fas] >                  # Fasta formatted genome file
#                                             < --offset [150] >                                                  # Minimum library size
#                                             < --distance [100] >                                             # Maximum variation above offset
#                                             [ --output [-] ]                                                          # Filename to write results to (default is STDOUT)
#                                             [ --readsize [35] ]                                                   # Raw sequence length
#                                             [ --usage ]                                                                # Prints this  
#
# Takes a paired ends sequencing output (needs to be preprocessed to fit above)
# and compares each pair to maximize unique alignment matches to a genome.
# Considers the following occurrences: 
# (1) both sides of pair have 0 matches; (2) one side of the pair has unique match;
# (3) one side of the pair has multiple matches; (4) both sides of the pair have unique matches;
# (5) one side of the pair has unique match, other has multiple matches; (6) both sides of the pair have multiple matches.
# Considers the following criteria:
# (a) Each side of a pair matches to a reverse complement chromosome of its own;
# (b) The distance between potential pairs in the case of multiple matches must be larger than $offset and smaller than $offset+$distance
#
# Last edited 2008-11-06
# Copyright 2008 Pedro Silva <psilva@dzlab.pmb.berkeley.edu/>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use diagnostics;
use Getopt::Long;

# Globals, passed as command line options
my $leftendfile   = ""; # left end sequences file in eland 3 format
my $rightendfile = ""; # right end sequences file in eland 3 format
my $referencefile = ""; # unmodified reference genome file in fasta format
my $distance = 100; # maximum (random) range distance beyond offset.
my $offset = 150; # minimum distance between left and right sequences of a pair.
my $output = '-'; # output mode. standard output by default
my $readsize = 35;
# For any given in silico sequence, the true distance was given by: (int(rand($distance)) + $offset + 1) --pedro
my $usage = 0; # print usage and exit

# Initial check of command line parameters
if ($usage || @ARGV<5) {
    print STDERR "correlatePairedEnds.pl\t<--left filename [leftendfile.eland3]>\t\t5' raw sequence file
\t\t\t<--right filename [rightendfile.eland3]>\t3' raw sequence file
\t\t\t<--reference [referencefile.fas]>\t\tFasta formatted genome file
\t\t\t<--offset [150]>\t\t\t\tMinimum library size
\t\t\t<--distance [100]>\t\t\t\tMaximum variation above offset
\t\t\t[--output [-]]\t\t\t\t\tFilename to write results to (default is STDOUT)
\t\t\t[--readsize [35]]\t\t\t\tRaw sequence length
\t\t\t[--usage]\t\t\t\t\tPrints this\nTakes a paired ends sequencing output (needs to be preprocessed to fit above) and compares each pair to maximize unique alignment matches to a genome.\n";
    exit 1;
}

# Parse command line arguments
my $result = GetOptions (
    "left|l=s" => \$leftendfile,
    "right|r=s" => \$rightendfile,
    "reference|ref=s"   => \$referencefile,
    "distance|d=i" => \$distance,
    "offset|off=i" => \$offset,
    "output|out:s" => \$output,
    "readsize|s:i" => \$readsize,
    "usage|help|h" => \$usage
    );

# holds name of chromosomes as keys and length of chromosomes in bp as values
my %reference=();

# begins limited scope for extracting names and lengths of chromosomes in reference file
{
    open(REF, "<$referencefile");
    my @fastaseq=<REF>;
    close(REF);

    my (@idx, @dsc)=();
    for(my $i=0;$i<@fastaseq;$i++) {
	if($fastaseq[$i] =~ m/^>/) {
	    $fastaseq[$i] =~ s/[>\n]//g;
	    $fastaseq[$i] = (split '\s', "$fastaseq[$i]")[0];
	    $fastaseq[$i] =~ tr/A-Z/a-z/;
	    push @idx, $i;
	    push @dsc, $fastaseq[$i];
	}
    }

    my $line;
    for(my $j=0;$j<@idx;$j++) {
	if($j==scalar(@idx)-1) {$line = join("", @fastaseq[$idx[$j]+1..@fastaseq-1]);}
	else {$line = join("", @fastaseq[$idx[$j]+1..$idx[$j+1]-1]);}
	$line=~s/\n//g;
	$reference{ $dsc[$j] } = length($line);
	$reference{ "$dsc[$j]-seq" } = $line;
	$reference{ "$dsc[$j]-rc" } = &reverseComp($line);
    }
}


# redirects STDOUT to file if specified by user
if(!($output eq '-')) {
    open(STDOUT, ">$output");
}
# reads left end sequences
open(LEFT, "<$leftendfile");
my @leftend=<LEFT>;
close(LEFT);
# reads right end sequences
open(RIGHT, "<$rightendfile");
my @rightend=<RIGHT>;
close(RIGHT);
# opens file for debugging
open(DBG, ">dbg");


# goes through left sequences
for(my $i=0;$i<@leftend;$i++) {

    my %left=%{&parseEland3Line($leftend[$i])};
    my %right=%{&parseEland3Line($rightend[$i])};

    my $lmatch=$left{'matches'};
    my $rmatch=$right{'matches'};

    my $lsequence=$left{'sequence'};
    my $rsequence=$right{'sequence'};

    my $allowed;
    $allowed=$left{'allowed'} or $allowed=$right{'allowed'} or $allowed=0;
    

    ##### NM - NM #####
    ##### No matches on either end #####
    if($lmatch==0 && $rmatch==0) {
	print $left{'line'} . "\tNM\t" . $lsequence . "\n";
	print $right{'line'} + @leftend . "\tNM\t" . $rsequence . "\n";
    }
    ##### No matches on either end #####



    ##### U0 - U0) #####
    ##### Unique matches on both ends #####
    if($lmatch==1 && $rmatch==1) {
	$left{'chr0'}=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	my $tmp=$1;
	$tmp=~tr/A-Z/a-z/;
	my $match=&checkMatch($left{'chr0'}, $right{'chr0'}, $left{'coord'}, $right{'coord'}, $reference{$tmp}, $offset, $distance);

	if($match) {
	    if($left{'chr0'}=~/^RC_/i) {
		print $left{'line'} . "\t" . "U$allowed" . "\t" . $lsequence . "\t" . ($reference{$tmp}-$left{'coord'}) . "\t" . substr($reference{"$tmp-rc"}, $left{'coord'}-1, $readsize) . "\t"  . abs($left{'coord'}-($reference{$tmp}-$right{'coord'}))  . "\n";
		print $right{'line'} + @rightend . "\t" . "U$allowed" .  "\t" . $rsequence . "\t" . $right{'coord'} . "\t" .  substr($reference{"$tmp-seq"}, $right{'coord'}-1, $readsize) . "\t"  . abs($left{'coord'}-($reference{$tmp}-$right{'coord'})) . "\n";
	    }
	    else {
		print $left{'line'} . "\t" .  "U$allowed" . "\t" . $lsequence . "\t" . $left{'coord'} . "\t" . substr($reference{"$tmp-seq"}, $left{'coord'}-1, $readsize) . "\t" . abs($left{'coord'}-($reference{$tmp}-$right{'coord'}))  ."\n";
		print $right{'line'} + @rightend . "\t" .  "U$allowed" . "\t" . $rsequence . "\t"  . ($reference{$tmp}-$right{'coord'}) . "\t" . substr($reference{"$tmp-rc"}, $right{'coord'}-1, $readsize) . "\t" . abs($left{'coord'}-($reference{$tmp}-$right{'coord'}))  . "\n";
	    }
	}
	else {
	    print $left{'line'} . "\tNM\t" . $lsequence . "\n";
	    print $right{'line'} + @rightend . "\tNM\t" . $rsequence . "\n";
	}
    }
    ##### Unique matches on both ends #####



    ##### NM - U0 #####
    ##### No match on one end, 1 match on other #####
    if($lmatch==0 && $rmatch==1) {
	$right{'chr0'}=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	my $tmp=$1;
	$tmp=~tr/A-Z/a-z/;
	if($right{'chr0'}=~m/^RC_/i) { # checks that this maps to reverse strand
	    my $revcoord=$reference{$tmp}-$right{'coord'};
	    $tmp=substr($reference{"$tmp-rc"}, $revcoord, $readsize);
	}
	else {
	    $tmp=substr($reference{"$tmp-seq"}, $right{'coord'}, $readsize);
	}
	print $left{'line'} . "\tNM\t" . $lsequence . "\n";
	print $right{'line'} + @rightend . "\t" . "U$allowed" . "\t" .  $rsequence . "\t" . ($reference{$tmp}-$right{'coord'}) . "\t" . $tmp . "\n";
    }

    if($lmatch==1 && $rmatch==0) {
	$left{'chr0'}=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	my $tmp=$1;
	$tmp=~tr/A-Z/a-z/;
	if($left{'chr0'}=~m/^RC_/i) {
	    my $revcoord=$reference{$tmp}-$left{'coord'};
	    $tmp=substr($reference{"$tmp-rc"}, $revcoord, $readsize);
	}
	else {
	    $tmp=substr($reference{"$tmp-seq"}, $left{'coord'}, $readsize);
	}
	print $left{'line'} + @leftend . "\t" . "U$allowed" . "\t" . $lsequence . "\t" .  ($reference{$tmp}-$left{'coord'}) . "\t" .  $tmp . "\n";
	print $right{'line'} . "\tNM\t" . $rsequence . "\n";
    }
    ##### No match on one end, 1 match on other #####



    ##### NM - R0 #####
    ##### No match on one end, multiple matches on other #####
    if($lmatch==0 && $rmatch>1) {
	print $left{'line'} . "\tNM\t" . $lsequence . "\n";
	print $right{'line'} + @rightend . "\t" . "R$allowed" . "\t" . $rsequence . "\t";
	for(my $i=0;$i<$rmatch;$i++) {
	    $right{"chr$i"}=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	    my $tmp=$1;
	    $tmp=~tr/A-Z/a-z/;
	    if($right{"chr$i"}=~m/^RC_/i) {
		print substr($reference{"$tmp-rc"}, $right{"coord$i"}-1, $readsize) . ",";
	    }
	    else {
		print substr($reference{"$tmp-seq"}, $right{"coord$i"}-1, $readsize) . ",";
	    }
	}
	print "\n";
    }

    if($lmatch>1 && $rmatch==0) {
	print $left{'line'} . "\t" . "R$allowed" . "\t" . $lsequence . "\t";
	for(my $i=0;$i<$lmatch;$i++) {
	    $left{"chr$i"}=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	    my $tmp=$1;
	    $tmp=~tr/A-Z/a-z/;
	    if($left{"chr$i"}=~m/^RC_/i) {
		print substr($reference{"$tmp-rc"}, $left{"coord$i"}-1, $readsize) . ",";
	    }
	    else {
		print substr($reference{"$tmp-seq"}, $left{"coord$i"}-1, $readsize) . ",";
	    }
	}
	print "\n";
	print $right{'line'} + @rightend . "\tNM\t" . $rsequence . "\n";
    }
    ##### No match on one end, multiple matches on other #####



    ##### One match on one end, multiple matches on other end #####
    if($lmatch==1 && $rmatch>1) {
	$left{'chr0'}=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	my $tmp=$1;
	$tmp=~tr/A-Z/a-z/;

 	my ($bestcoord, $beststrand);
	my $score=100000000;
 	for(my $i=0;$i<$rmatch;$i++) {
	    my $tmpscore=&checkMultMatch($left{'chr0'}, $right{"chr$i"}, $left{'coord'}, $right{"coord$i"}, $reference{$tmp}, $offset, $distance);
	    if($tmpscore==-1) {next;}
	    if($tmpscore < $score) {
		$score=$tmpscore;
		$bestcoord=$right{"coord$i"};
		$beststrand=$right{"chr$i"};
	    }
	    else {next;}
 	}

	$beststrand=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	$tmp=$1;
	$tmp=~tr/A-Z/a-z/;
	if($beststrand=~m/^RC_/i) {
	    print $left{'line'} . "\t" . "U$allowed" . "\t" . $lsequence . "\t" . $left{'coord'} . "\t" .  substr($reference{"$tmp-seq"}, $left{'coord'}-1, $readsize) . "\t" . $score . "\n";
	    print $right{'line'} + @rightend . "\t" . "U$allowed" . "\t" . $rsequence . "\t" .  ($reference{$tmp}-$bestcoord) . "\t" .  substr($reference{"$tmp-rc"}, $bestcoord-1, $readsize) . "\t" . $score .  "\n";
	}
	else {
	    print $left{'line'} . "\t" . "U$allowed" . "\t" . $lsequence . "\t" .  ($reference{$tmp}-$left{'coord'}) . "\t" . substr($reference{"$tmp-rc"}, $left{'coord'}-1, $readsize) . "\t" . $score .  "\n";
	    print $right{'line'} + @rightend . "\t" . "U$allowed" . "\t" . $rsequence . "\t" . $bestcoord . "\t" . substr($reference{"$tmp-seq"}, $bestcoord-1, $readsize) . "\t" . $score .  "\n";
	}
    }

    if($lmatch>1 && $rmatch==1) {
	$right{'chr0'}=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	my $tmp=$1;
	$tmp=~tr/A-Z/a-z/;

 	my ($bestcoord, $beststrand);
	my $score=100000000;
 	for(my $i=0;$i<$lmatch;$i++) {
	    my $tmpscore=&checkMultMatch($left{"chr$i"}, $right{"chr0"}, $left{"coord$i"}, $right{"coord"}, $reference{$tmp}, $offset, $distance);
	    if($tmpscore==-1) {next;}
	    if($tmpscore < $score) {
		$score=$tmpscore;
		$bestcoord=$left{"coord$i"};
		$beststrand=$left{"chr$i"};
	    }
 	}

	$beststrand=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	$tmp=$1;
	$tmp=~tr/A-Z/a-z/;
	if($beststrand=~m/^RC_/i) {
	    print $left{'line'} . "\t" . "U$allowed" . "\t" . $lsequence . "\t" .  ($reference{$tmp}-$bestcoord) . "\t" . substr($reference{"$tmp-rc"}, $bestcoord-1, $readsize) . "\t"  . $score .  "\n";
	    print $right{'line'} + @rightend . "\t" . "U$allowed" . "\t" . $rsequence . "\t" . $right{'coord'} . "\t" . substr($reference{"$tmp-seq"}, $right{'coord'}-1, $readsize) . "\t" . $score .  "\n";
	}
	else {
	    print $left{'line'} . "\t" . "U$allowed" . "\t" . $lsequence . "\t" . $bestcoord . "\t" . substr($reference{"$tmp-seq"}, $bestcoord-1, $readsize) . "\t" . $score .  "\n";
	    print $right{'line'} + @rightend . "\t" . "U$allowed" . "\t" . $rsequence . "\t" .  ($reference{$tmp}-$right{'coord'}) . "\t" . substr($reference{"$tmp-rc"}, $right{'coord'}-1, $readsize) . "\t" . $score .  "\n";
	}
    }
   #####    One match on one end, multiple matches on other end #####



    ##### R0 - R0 #####
    ##### Multiple matches on both ends #####
    if($lmatch>1 && $rmatch>1) {

	print DBG "\n==============================\n";
	print DBG "SEQ ", $i+1, "\n";

	my ($lbestcoord, $rbestcoord, $lbeststrand, $rbeststrand);
	my $bestscore=100000000;

	# loops through each possible target on the left
	for(my $i=0;$i<$lmatch;$i++) {

	    $left{"chr$i"}=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	    my $tmp=$1;
	    $tmp=~tr/A-Z/a-z/;

	    my ($bestcoord, $beststrand);
	    my $score=100000000;

	    print DBG "lmatch: ", $i, "\n";

	    # loops through each possible target on the right
	    for(my $j=0;$j<$rmatch;$j++) {

		my $tmpscore=&checkMultMatch($left{"chr$i"}, $right{"chr$j"}, $left{"coord$i"}, $right{"coord$j"}, $reference{$tmp}, $offset, $distance);

		if($tmpscore==-1) {next;}
		if($tmpscore < $score) {
		    $score=$tmpscore;
		    $bestcoord=$right{"coord$j"};
		    $beststrand=$right{"chr$j"};
		}
		print DBG "\trmatch: ", $j, "\tright chr: ", $right{"chr$j"}, "\tbestcoord: ", $bestcoord, "\tbeststrand: ", $beststrand, "\tscore: ", $score, "\n";
	    }
 
	    print DBG "\nltarget: ",  $left{"chr$i"}, "\tbestcoord: ", $bestcoord, "\tbeststrand: ", $beststrand, "\tscore: ", $score, "\n\n";

	    # if no good matches on right side, skip this left sequence
	    if($bestcoord==-1 && $beststrand==-1) {next;}
 
	    if($score < $bestscore) {
		$bestscore=$score;
		$lbestcoord=$left{"coord$i"};
		$lbeststrand=$left{"chr$i"};
		$rbestcoord=$bestcoord;
		$rbeststrand=$beststrand;
	    }
	}
	
	$lbeststrand=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	my $tmp=$1;
	$tmp=~tr/A-Z/a-z/;
	if($lbeststrand=~m/^RC_/i) {
	    print $left{'line'} . "\t" . "U$allowed" . "\t" .  $lsequence . "\t" .  ($reference{$tmp}-$lbestcoord) . "\t" . substr($reference{"$tmp-rc"}, $lbestcoord-1, $readsize) . "\t" . $bestscore . "\n";
	    print $right{'line'} + @rightend . "\t" . "U$allowed" . "\t" .  $rsequence . "\t" . $rbestcoord . "\t" . substr($reference{"$tmp-seq"}, $rbestcoord-1, $readsize) . "\t" . $bestscore .  "\n";
	}
	else {
	    print $left{'line'} . "\t" . "U$allowed" . "\t" . $lsequence . "\t" . $lbestcoord . " \t" . substr($reference{"$tmp-seq"}, $lbestcoord-1, $readsize) . "\t" . $bestscore .  "\n";
	    print $right{'line'} + @rightend . "\t"  . "U$allowed" . "\t" . $rsequence . "\t" .  ($reference{$tmp}-$rbestcoord) . "\t" .  substr($reference{"$tmp-rc"}, $rbestcoord-1, $readsize) . "\t" . $bestscore .  "\n";
	}
    }
    ##### Multiple matches on both ends #####
}
# end for loop through all sequences
close(DBG);
# main program finished


##### Start function definitions #####

# Given 2 chromosomes and 2 coordinates, plus a library offset and distance (see header for details) checks
# that each element of a pair belongs to the reverse complement of the other, and that both their coordinates
# summed don't exceed the length of the chromosome plus or minus the library size.
# Returns $match=1  or 0
sub checkMatch {
    my ($lchr, $rchr, $lcoord, $rcoord, $chrlength, $offset, $distance)=@_;
    my $match=1;

    # check chromosome
    if($lchr=~m/^(chr.)/i) {
	if(!$rchr=~m/^RC_$1/i) {$match=0;}
    }
    else {
	if($rchr=~m/^RC_(chr.)/i) {
	    if(!$lchr=~m/^$1/i) {$match=0;}
	}
    }

    # check coordinates
    if( ($lcoord+$rcoord>$chrlength+($offset+$distance)) || ($lcoord+$rcoord<$chrlength-($offset+$distance)) ) {
	$match=0;
    }
    return $match;
}

# Given 2 chromosomes and 2 coordinates, plus a library offset and distance (see header for details) checks
# that each element of a pair belongs to the reverse complement of the other, and that both their coordinates
# summed don't exceed the length of the chromosome plus or minus the library size.
# Returns $score=library size/distance between both sequences
sub checkMultMatch {
    my ($lchr, $rchr, $lcoord, $rcoord, $chrlength, $offset, $distance)=@_;
    my $score;
    
    # check chromosome
    if($lchr=~m/^(chr.)/i) {
	if(!$rchr=~m/^RC_$1/i) {$score=-1;}
    }
    else {
	if($rchr=~m/^RC_(chr.)/i) {
	    if(!$lchr=~m/^$1/i) {$score=-1;}
	}
    }
    # check absolute coordinates
    if( ($lcoord+$rcoord>$chrlength+($offset+$distance)) || ($lcoord+$rcoord<$chrlength-($offset+$distance)) ) {
	$score=-1;
    }
    # check relative coordinates    
    $score=abs($lcoord-($chrlength-$rcoord));
    return $score;
}

# Parses single eland3 record (output from seqmap)
# Returns reference to hash with following keys:
# 'line', 'sequence', 'matches', 'chr[0-#matches-1]', 'coord[0-#matches-1]'
sub parseEland3Line {
    my (@line)=(split(/\t/, shift));

    my %hash=();
    $hash{'line'}=$line[0];
    $hash{'sequence'}=$line[1];
    if($line[2] =~ 'NM') {$hash{'matches'}=0;}
    else {$hash{'matches'}=(split ':', $line[2])[0];}
    my @all=split(/,/, $line[3]);
    $hash{'matches'}=scalar(@all);
    if($hash{'matches'}==1) {
	my $tmp=(split ':', $line[3])[0];
	$tmp=~s/\n//g;
	$hash{'chr0'}=$tmp;
	$tmp=(split ':', $line[3])[1];
	$tmp=~s/[A-Z]([0-9])$//i;
	$hash{'allowed'}=$1;
	$tmp=~s/\n//g;
	$hash{'coord'}=$tmp;
    }
    if($hash{'matches'}>1) {
	my @all=split(/,/, $line[3]);
	$hash{'matches'}=scalar(@all);
	for(my $i=0;$i<@all;$i++) {
	    my $tmp=(split ':', $all[$i])[0];
	    $tmp=~s/\n//g;
	    $hash{"chr$i"}=$tmp;
	    $tmp=(split ':', $all[$i])[1];
	    $tmp=~s/[A-Z][0-9]$//i;
	    $tmp=~s/\n//g;
	    $hash{"coord$i"}=$tmp;
	}
    }
    return \%hash;
}

# Converts input sequence to reverse complement
# Returns scalar $string with processed sequence
sub reverseComp {
    my $shortseq = $_[0];
    $shortseq =~ tr/ACGTacgt/TGCAtgca/;
    $shortseq =~ s/\n//;
    return reverse $shortseq;
}
