#!/usr/bin/perl
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
#

=head1 Synopsis

Takes a paired ends sequencing output (needs to be preprocessed to fit above)
and compares each pair to maximize unique alignment matches to a genome.
Considers the following occurrences:
(1) both sides of pair have 0 matches; (2) one side of the pair has unique match;
(3) one side of the pair has multiple matches; (4) both sides of the pair have unique matches;
(5) one side of the pair has unique match, other has multiple matches; (6) both sides of the pair have multiple matches.
Considers the following criteria:
(a) Each side of a pair matches to a reverse complement chromosome of its own;
(b) The distance between potential pairs in the case of multiple matches must be larger than $offset and smaller than $offset+$distance
=cut

use strict;
use warnings;
#use diagnostics;
use Getopt::Long;
use threads;

# Globals, passed as command line options
my $leftendfile   = ""; # left end sequences file in eland 3 format
my $rightendfile = ""; # right end sequences file in eland 3 format
my $referencefile = ""; # unmodified reference genome file in fasta format
my $distance = 150; # maximum (random) range distance beyond offset.
my $offset = 50; # minimum distance between left and right sequences of a pair.
my $output = '-'; # output mode. standard output by default
my $readsize = 45;
my $output_repeats = 0;
# For any given in silico sequence, the true distance was given by: (int(rand($distance)) + $offset + 1) --pedro
my $usage = 0; # print usage and exit

# Initial check of command line parameters
if ($usage || @ARGV<5) {
    print STDERR
"correlatePairedEnds.pl <PARAMETERS> [OPTIONS]
\t<--left>\t5' raw sequence file
\t<--right>\t3' raw sequence file
\t<--reference>\tFasta genome file
\t<--offset>\tMinimum library size
\t<--distance>\tMaximum variation from offset
\t[--output]\tFilename to write results to (default is STDOUT)
\t[--readsize]\tRaw sequence length
\t[--repeats]\tPrint reads with multiple matches
\t[--usage]\tPrints this
Takes a paired ends sequencing output (needs to be preprocessed to fit above) 
and compares each pair to maximize unique alignment matches to a genome.\n";
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
    "repeats" => \$output_repeats,
    "usage|help|h" => \$usage
    );

# holds name of chromosomes as keys and length of chromosomes in bp as values
my %reference=();

# begins limited scope for extracting names and lengths of chromosomes in reference file
{
    print STDERR ("Indexing reference file...");
    open(my $REF, "<", "$referencefile");
    my @fastaseq=<$REF>;
    close($REF);

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
    print STDERR "OK.";
}


# redirects STDOUT to file if specified by user
if(!($output eq '-')) {
    open(STDOUT, ">", "$output");
}

# reads left end sequences
open(my $LEFT, "<", "$leftendfile");
my @leftend=<$LEFT>;
close($LEFT);

# reads right end sequences
open(my $RIGHT, "<", "$rightendfile");
my @rightend=<$RIGHT>;
close($RIGHT);

print STDERR "\nCorrelating pairs";

# goes through left sequences
for(my $i=0;$i<@leftend;$i++) {

    print STDERR "." if(int(rand(@leftend)<@leftend/1000));

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
	print join("\t",
		   $left{'line'},
		   "NM",
		   $lsequence,
		   "\n");
	print join("\t",
		   $right{'line'} + @leftend,
		   "NM",
		   $rsequence,
		   "\n");
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
		print join("\t",
			   $left{'line'},
			   "U$allowed",
			   $lsequence,
			   ($reference{$tmp}-$left{'coord'}),
			   substr($reference{"$tmp-rc"}, $left{'coord'}-1, $readsize+2),
			   abs($left{'coord'}-($reference{$tmp}-$right{'coord'})),
			   "\n");
		print join("\t",
			   $right{'line'} + @rightend,
			   "U$allowed",
			   $rsequence,
			   $right{'coord'},
			   substr($reference{"$tmp-seq"}, $right{'coord'}-1, $readsize+2),
			   abs($left{'coord'}-($reference{$tmp}-$right{'coord'})),
			   "\n");
	    }
	    else {
		print join("\t",
			    $left{'line'},
			    "U$allowed",
			    $lsequence,
			    $left{'coord'},
			    substr($reference{"$tmp-seq"}, $left{'coord'}-1, $readsize+2),
			    abs($left{'coord'}-($reference{$tmp}-$right{'coord'})),
			    "\n");
		print join("\t",
			    $right{'line'} + @rightend,
			    "U$allowed",
			    $rsequence,
			    ($reference{$tmp}-$right{'coord'}),
			    substr($reference{"$tmp-rc"}, $right{'coord'}-1, $readsize+2),
			    abs($left{'coord'}-($reference{$tmp}-$right{'coord'})),
			    "\n");
	    }
	}
	else {
	    print join("\t",
			$left{'line'},
			"NM",
			$lsequence,
			"\n");
	    print join("\t",
		       $right{'line'} + @rightend,
		       "NM",
		       $rsequence,
		       "\n");
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
	    $tmp=substr($reference{"$tmp-rc"}, $revcoord, $readsize+2);
	}
	else {
	    $tmp=substr($reference{"$tmp-seq"}, $right{'coord'}, $readsize+2);
	}
	print join("\t",
		   $left{'line'},
		   "NM",
		   $lsequence,
		   "\n");
	print join("\t",
		   $right{'line'} + @rightend,
		   "U$allowed",
		   $rsequence,
		   ($reference{$tmp}-$right{'coord'}),
		   $tmp,
		   "\n");
    }

    if($lmatch==1 && $rmatch==0) {
	$left{'chr0'}=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	my $tmp=$1;
	$tmp=~tr/A-Z/a-z/;
	if($left{'chr0'}=~m/^RC_/i) {
	    my $revcoord=$reference{$tmp}-$left{'coord'};
	    $tmp=substr($reference{"$tmp-rc"}, $revcoord, $readsize+2);
	}
	else {
	    $tmp=substr($reference{"$tmp-seq"}, $left{'coord'}, $readsize+2);
	}
	print join("\t",
		   $left{'line'} + @leftend,
		   "U$allowed",
		   $lsequence,
		   ($reference{$tmp}-$left{'coord'}),
		   $tmp,
		   "\n");
	print join("\t",
		   $right{'line'},
		   "NM",
		   $rsequence,
		   "\n");
    }
    ##### No match on one end, 1 match on other #####



    ##### NM - R0 #####
    ##### No match on one end, multiple matches on other #####
    if($output_repeats) {
	if($lmatch==0 && $rmatch>1) {
	    print join("\t",
		       $left{'line'},
		       "NM",
		       $lsequence,
		       "\n");
	    print join("\t",
		       $right{'line'} + @rightend,
		       "R$allowed",
		       $rsequence,
		       "\t");
	    for(my $i=0;$i<$rmatch;$i++) {
		$right{"chr$i"}=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
		my $tmp=$1;
		$tmp=~tr/A-Z/a-z/;
		if($right{"chr$i"}=~m/^RC_/i) {
		    print substr($reference{"$tmp-rc"}, $right{"coord$i"}-1, $readsize+2) . ",";
		}
		else {
		    print substr($reference{"$tmp-seq"}, $right{"coord$i"}-1, $readsize+2) . ",";
		}
	    }
	    print "\n";
	}

	if($lmatch>1 && $rmatch==0) {
	    print join("\t",
		       $left{'line'},
		       "R$allowed",
		       $lsequence,
		       "\t");
	    for(my $i=0;$i<$lmatch;$i++) {
		$left{"chr$i"}=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
		my $tmp=$1;
		$tmp=~tr/A-Z/a-z/;
		if($left{"chr$i"}=~m/^RC_/i) {
		    print substr($reference{"$tmp-rc"}, $left{"coord$i"}-1, $readsize+2) . ",";
		}
		else {
		    print substr($reference{"$tmp-seq"}, $left{"coord$i"}-1, $readsize+2) . ",";
		}
	    }
	    print "\n";
	    print join("\t",
		       $right{'line'} + @rightend,
		       "NM",
		       $rsequence,
		       "\n");
	}
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
	    print join("\t",
		       $left{'line'},
		       "U$allowed",
		       $lsequence,
		       $left{'coord'},
		       substr($reference{"$tmp-seq"}, $left{'coord'}-1, $readsize+2),
		       $score,
		       "\n");
	    print join("\t",
		       $right{'line'} + @rightend,
		       "U$allowed",
		       $rsequence,
		       ($reference{$tmp}-$bestcoord),
		       substr($reference{"$tmp-rc"}, $bestcoord-1, $readsize+2),
		       $score,
		       "\n");
	}
	else {
	    print join("\t",
		       $left{'line'},
		       "U$allowed",
		       $lsequence,
		       ($reference{$tmp}-$left{'coord'}),
		       substr($reference{"$tmp-rc"}, $left{'coord'}-1, $readsize+2),
		       $score,
		       "\n");
	    print join("\t",
		       $right{'line'} + @rightend,
		       "U$allowed",
		       $rsequence,
		       $bestcoord,
		       substr($reference{"$tmp-seq"}, $bestcoord-1, $readsize+2),
		       $score,
		       "\n");
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
	    print join("\t",
		       $left{'line'},
		       "U$allowed",
		       $lsequence,
		       ($reference{$tmp}-$bestcoord),
		       substr($reference{"$tmp-rc"}, $bestcoord-1, $readsize+2),
		       $score,
		       "\n");
	    print join("\t",
		       $right{'line'} + @rightend,
		       "U$allowed",
		       $rsequence,
		       $right{'coord'},
		       substr($reference{"$tmp-seq"}, $right{'coord'}-1, $readsize+2),
		       $score,
		       "\n");
	}
	else {
	    print join("\t",
		       $left{'line'},
		       "U$allowed",
		       $lsequence,
		       $bestcoord,
		       substr($reference{"$tmp-seq"}, $bestcoord-1, $readsize+2),
		       $score,
		       "\n");
	    print join("\t",
		       $right{'line'} + @rightend,
		       "U$allowed",
		       $rsequence,
		       ($reference{$tmp}-$right{'coord'}),
		       substr($reference{"$tmp-rc"}, $right{'coord'}-1, $readsize+2),
		       $score,
		       "\n");
	}
    }
   #####    One match on one end, multiple matches on other end #####



    ##### R0 - R0 #####
    ##### Multiple matches on both ends #####
    if($lmatch>1 && $rmatch>1) {
	my ($lbestcoord, $rbestcoord, $lbeststrand, $rbeststrand);
	my $bestscore=100000000;

	# loops through each possible target on the left
	for(my $i=0;$i<$lmatch;$i++) {

	    $left{"chr$i"}=~m/(chr.)/i; # this puts the base chromosome name ('chr.') into $1
	    my $tmp=$1;
	    $tmp=~tr/A-Z/a-z/;

	    my ($bestcoord, $beststrand);
	    my $score=100000000;

	    # loops through each possible target on the right
	    for(my $j=0;$j<$rmatch;$j++) {

		my $tmpscore=&checkMultMatch($left{"chr$i"}, $right{"chr$j"}, $left{"coord$i"}, $right{"coord$j"}, $reference{$tmp}, $offset, $distance);

		if($tmpscore==-1) {next;}
		if($tmpscore < $score) {
		    $score=$tmpscore;
		    $bestcoord=$right{"coord$j"};
		    $beststrand=$right{"chr$j"};
		}
	    }

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
	    print join("\t",
		       $left{'line'},
		       "U$allowed",
		       $lsequence,
		       ($reference{$tmp}-$lbestcoord),
		       substr($reference{"$tmp-rc"}, $lbestcoord-1, $readsize+2),
		       $bestscore,
		       "\n");
	    print join("\t",
		       $right{'line'} + @rightend,
		       "U$allowed",
		       $rsequence,
		       $rbestcoord,
		       substr($reference{"$tmp-seq"}, $rbestcoord-1, $readsize+2),
		       $bestscore,
		       "\n");
	}
	else {
	    print join("\t",
		       $left{'line'},
		       "U$allowed",
		       $lsequence,
		       $lbestcoord,
		       substr($reference{"$tmp-seq"}, $lbestcoord-1, $readsize+2),
		       $bestscore,
		       "\n");
	    print join("\t",
		       $right{'line'} + @rightend,
		       "U$allowed",
		       $rsequence,
		       ($reference{$tmp}-$rbestcoord),
		       substr($reference{"$tmp-rc"}, $rbestcoord-1, $readsize+2),
		       $bestscore,
		       "\n");
	}
    }
    ##### Multiple matches on both ends #####
}
# end for loop through all sequences
print STDERR "done!\n";
close(STDOUT);

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
    if( ($lcoord+$rcoord>$chrlength+($offset+$distance)*2) || ($lcoord+$rcoord<$chrlength-($offset+$distance)*2) ) {
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
    if( ($lcoord+$rcoord>$chrlength+($offset+$distance)*2) || ($lcoord+$rcoord<$chrlength-($offset+$distance)*2) ) {
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

    if($line[2] =~ 'NM') {
	$hash{'matches'}=0;
    }
    else {
	$hash{'matches'}=(split ':', $line[2])[0];
    }

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
