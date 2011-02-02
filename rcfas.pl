#!/usr/bin/env perl
# ___UNDOCUMENTED___
#
# rcfas.pl <inputfile> [> outputfile]
#       where: <inputfile> is a fasta formated sequence file
# Reads a sequence or set of sequences in fasta format
# Cleans up sequence header (ie, for each chromosome)
# Computes reverse complement for each sequence
# Outputs to standard output, so use redirection to save.
#
# Last edited 2008-10-14
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

unless (@ARGV == 1) {
    print STDERR ("Usage: rc.pl <fastafile.fas>\n") ;
    exit(-1);
}

# reads fasta file contents to array
my $fastafile = "$ARGV[0]";
open(FAS,"<$fastafile") or die "Can't read input file";
my @fastaseq=<FAS>;
close(FAS);

# @idx is array with indices for starting line of each sequence
# @dsc is array with the comment line ('>*') for each sequence
# @chr is the actual sequence content
my (@idx, @dsc, @chr)=();

# goes through each line in input array
# cleans up new lines, grab only first identifier and add RC prefix
# saves line number of each new chromosome in @idx and contents in @dsc
for(my $i=0;$i<@fastaseq;$i++) {
    if($fastaseq[$i] =~ m/^>/) {
	$fastaseq[$i] =~ s/[\n\r]//g;
	$fastaseq[$i] = (split '\s', "$fastaseq[$i]")[0];
	$fastaseq[$i] .= "\n";
	my $line = $fastaseq[$i];
	$line =~ s/^>(.+)/>RC_$1/i;
	push @idx, $i;
	push @dsc, $line;
    }
}	

# goes through every chromosome
# puts each chromosome into $line, deletes new lines, complements and reverses $line
# puts 80 characters of sequence into @chr at a time, after each description
for(my $j=0;$j<@idx;$j++) {
    my $line;
    if($j>0){push @chr, "\n";}
    push @chr, $dsc[$j];
    if($j==scalar(@idx)-1) {$line = join("", @fastaseq[$idx[$j]+1..@fastaseq-1]);}
    else {$line = join("", @fastaseq[$idx[$j]+1..$idx[$j+1]-1]);}
    $line =~ s/[\n\r]//g;
    $line =~ tr/ACGTacgt/TGCAtgca/;
    $line = reverse $line;
    my $numlines=length($line)/80;
    if((length($line) % 80) != 0) {$numlines++;}

    for(my $k=0;$k<$numlines;$k++) {
	if($k*80>length($line)) {next;}
	if($k*80+80>length($line)) {
	    my $tmp = substr($line, ($k * 80));
	    push @chr, $tmp;
	}
	else {
	    my $tmp = substr($line, ($k * 80), 80);
	    $tmp .= "\n";
	    push @chr, $tmp;
	}
    }
}

# prints original and RC sequences to standard output
print @fastaseq;
print @chr;

exit(0);
