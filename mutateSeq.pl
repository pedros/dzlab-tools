#!/usr/bin/env perl
#
# mutateSeq.pl <inputfile> [numseqs] [seqsize]
#       where: <inputfile> is a fasta formated sequence file
#              [numseqs] is an the number of random short sequences to output
#              [seqsize] is the width in base pairs of each short sequence.
# Simulates bisulfite treatment of DNA by converting cytosine to thymine in 90% of all sequences, with 10% probability of converting remaining unchanged sequences
#
# Last edited 2008-10-02
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

use warnings;
use strict;
use List::Util qw/shuffle/;  

die "Usage:\tmutateSeq.pl <inputfile> [numseqs] [seqsize]\n\t<inputfile> is a fasta formated sequence file\n\t[numseqs] is an the number of random short sequences to output\n\t[seqsize] is the width in base pairs of each short sequence.\nSimulates bisulfite treatment of DNA by converting cytosine to thymine in 90% of all sequences, with 10% probability of converting remaining unchanged sequences.\n"
    unless(@ARGV>0);

##### START MAIN #####

# single or paired ends output file
my ($paired, $distance, $offset)=(1, 150, 50);

# read fasta file contents to string
my $fastafile = "$ARGV[0]";
my $fastaseq = &readFAS($fastafile);

# get $numseqs random seqs of $shortseqsize
my $numseqs=10000;
my $shortseqsize = 35;
$numseqs=$ARGV[1] if($ARGV[1]);
$shortseqsize=$ARGV[2] if($ARGV[2]);
my @shortseqs = &fas2short($numseqs, $shortseqsize, $fastaseq);

# compute reverse complement for random half the short sequences
# we can do it linearly because fas2short randomizes picks
my $numrevcomp = int($numseqs * 0.5);
for(my $i=0;$i<$numrevcomp;$i++) {
    # reverse complement for single ends output
    if(!$paired) {$shortseqs[$i]=&reverseComp($shortseqs[$i]) . "\n";}
    # reverse complement for paired ends output
    # reverses the whole pair
    else {
        my @pairedends=split(/\t/, $shortseqs[$i]);
	$pairedends[0]=&reverseComp($pairedends[0]);
	$pairedends[1]=&reverseComp($pairedends[1]) . "\n";
	$shortseqs[$i]=join("\t", @pairedends);
    }
}

# shuffle to make sure we mutate normal and reverse complement sequences randomly
@shortseqs = shuffle(@shortseqs);

# mutate 90% of the sequences and mutating remaining ones with probability = 0.1
my $numC2T=int($numseqs * 0.9);
for(my $i=0;$i<$numseqs;$i++) {
    if($i<$numC2T) {
	$shortseqs[$i]=&mutateSeq($shortseqs[$i]);
    }
    else {
	if (int(rand(9))==1) {$shortseqs[$i]=&mutateSeq($shortseqs[$i]);}
    }
}

# Second reverse complement for paired ends output
if($paired) {
    for(my $i=0;$i<$numseqs;$i++) {
	# reverse complement for single ends output
	if(!$paired) {$shortseqs[$i]=&reverseComp($shortseqs[$i]) . "\n";}
	# reverse complement for paired ends output
	# on second end sequences only
	else {
	    my @pairedends=split(/\t/, $shortseqs[$i]);
	    $pairedends[1]=&reverseComp($pairedends[1]) . "\n";
	    $shortseqs[$i]=join("\t", @pairedends);
	}
    }
}

# final randomize to make sure results are not sorted by probability of mutation
@shortseqs = shuffle(@shortseqs);

print @shortseqs;

exit(0);

##### END MAIN #####


##### START SUBS #####

# reads a fasta formated sequence,
# deletes sequence headers and white space, and
# concatenates output to string
sub readFAS {
    open(FAS,"<$_[0]") or die "Can't read input file";
    while(my $line = <FAS>) {
	next unless $line !~ m/^>/;
	$line =~ s/\s//;
	$fastaseq .= $line;
    }
    close(FAS);
    return $fastaseq;
}

# converts a raw string containing a sequence to
# an array of shorter sequences of length $shortseqsize in random order
# with overlaps (p=$shortseqsize/$numseqs) or repetitions (p=1/$numseqs) allowed
sub fas2short {
    my @shortseq = ();
    my ($numseqs, $shortseqsize, $fastaseq)=($_[0], $_[1], $_[2]);
    for(my $i=0;$i<$numseqs;$i++) {

	my $shortseqstart = int(rand(length($fastaseq)));
	my $shortseqend = $shortseqstart + $shortseqsize;

	# grab single ends sequence
	if($shortseqend < length($fastaseq)) {
	    $shortseq[$i] = substr($fastaseq, $shortseqstart, $shortseqsize);
	}
	else {$shortseq[$i] = substr($fastaseq, $shortseqstart);}
    
	# grab paired ends sequence
	if($paired) {
	    $shortseq[$i] .= "\t";
	    my $secseqstart = $shortseqstart + (int(rand($distance)) + $offset + 1);
	    my $secseqend = $secseqstart + $shortseqsize;
	    if($secseqend < length($fastaseq)) {
		$shortseq[$i] .= substr($fastaseq, $secseqstart, $shortseqsize);
	    }
	    else {$shortseq[$i] .= substr($fastaseq, $secseqstart);}
	}
	$shortseq[$i] .= "\n";
    }
    return @shortseq;
}

# converts input sequence to reverse complement
sub reverseComp {
    my $shortseq = $_[0];
    $shortseq =~ tr/ACGTacgt/TGCAtgca/;
    $shortseq =~ s/\n//;
    return reverse $shortseq;
}

# simulates bisulfite treatment on cytosine
sub mutateSeq {
    my $shortseq = $_[0];
    $shortseq =~ tr/Cc/Tt/; 
    return $shortseq;
}

##### END SUBS #####
