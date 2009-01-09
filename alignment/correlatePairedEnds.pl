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

use Getopt::Long;
use strict;
use warnings;
use diagnostics;disable diagnostics;
use Carp;

# Globals, passed as command line options
my $leftendfile; # left end sequences file in eland 3 format
my $rightendfile; # right end sequences file in eland 3 format
my $referencefile; # unmodified reference genome file in fasta format
my $distance = 300; # maximum (random) range distance beyond offset.
my $offset = 0; # minimum distance between left and right sequences of a pair.
my $output = q{-}; # output mode. standard output by default
my $readsize = 45; # length of each read
my $repeats = 0; # print ambiguous reads
my $nomatches = 0; # print unmatched reads
my $usage = 0; # print usage and exit

# Initial check of command line parameters
usage ();

# Parse command line arguments
my $result = GetOptions (
    "left|l=s" => \$leftendfile,
    "right|r=s" => \$rightendfile,
    "reference|ref=s"   => \$referencefile,
    "distance|d=i" => \$distance,
    "offset|off=i" => \$offset,
    "output|out:s" => \$output,
    "readsize|s:i" => \$readsize,
    "repeats" => sub {$repeats=1;},
    "nomatches" => sub {$nomatches=1},
    "verbose|v" => sub {enable diagnostics;},
    "quiet|q" => sub {no warnings;},
    "usage|help|h" => \&usage
    );

# holds name of chromosomes as keys and length of chromosomes in bp as values
my %reference=();

# begins limited scope for extracting names and lengths of chromosomes in reference file
{
    print STDERR 'Indexing reference file...';
    open my $REF, '<', "$referencefile" or croak "Can't open file: $referencefile";
    my @fastaseq=<$REF>;
    close($REF);

    my (@idx, @dsc)=();
    for(my $i=0;$i<@fastaseq;$i++) {
	if($fastaseq[$i] =~ m/^>/) {
	    $fastaseq[$i] =~ s/>//g;
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
	$line=~s/[\n\r]//g;
	$reference{ $dsc[$j] } = length($line);
	$reference{ "$dsc[$j]-seq" } = $line;
	$reference{ "$dsc[$j]-rc" } = &reverseComp($line);
    }
    print STDERR "OK";
}


# redirects STDOUT to file if specified by user
if(!($output eq '-')) {
    open(STDOUT, ">", "$output");
}

# reads left end sequences
open my $LEFT, '<', $leftendfile or croak "Can't open file: $leftendfile";
#my @leftend=<$LEFT>;
#close($LEFT);

# reads right end sequences
open my $RIGHT, '<', $rightendfile  or croak "Can't open file: $rightendfile";
#my @rightend=<$RIGHT>;
#close($RIGHT);

print STDERR "\nCorrelating pairs...";

# goes through left sequences
while (<$LEFT>) {

    my $leftend = <$LEFT>;
    my $rightend = <$RIGHT>;

    my %left=%{&parseEland3Line($leftend)};
    my %right=%{&parseEland3Line($rightend)};

# for(my $i=0;$i<@leftend;$i++) {

#     my %left=%{&parseEland3Line($leftend[$i])};
#     my %right=%{&parseEland3Line($rightend[$i])};

    my $lmatch = $left{'matches'};
    my $rmatch = $right{'matches'};

    my $lsequence = $left{'sequence'};
    my $rsequence = $right{'sequence'};

    $lsequence =~ s/[\r\n]//g;
    $rsequence =~ s/[\r\n]//g;

    my ($l_seqname, $l_source, $l_feature, $l_start, $l_end, $l_score, $l_strand, $l_frame, $l_attribute) = ('.', 'pair_align', '.', 0, 0, 0, '.', '.', '.');
    my ($r_seqname, $r_source, $r_feature, $r_start, $r_end, $r_score, $r_strand, $r_frame, $r_attribute) = ('.', 'pair_align', '.', 0, 0, 0, '.', '.', '.');

    ##### NM - NM #####
    ##### No matches on either end #####
    if ($nomatches) {
	if ($lmatch == 0 && $rmatch == 0) {

            $l_feature   = $left{'line'} . ':' . $lsequence;
            $r_feature   = $right{'line'} . ':' . $rsequence;

            ### to delete
	    print join("\t",
		       '.',
		       'pair_align',
		       $left{'line'} . ':' . $lsequence,
		       '.',
		       '.',
		       0,
		       '.',
		       ".",
		       ".") . "\n";
	    print join("\t",
		       ".",
		       "pair_align",
		       $right{'line'} . ":" . $rsequence,
		       ".",
		       ".",
		       0,
		       ".",
		       ".",
		       ".") . "\n";
            ### to delete
	}
    }
    ##### No matches on either end #####



    ##### U0 - U0) #####
    ##### Unique matches on both ends #####
    if ($lmatch == 1 && $rmatch == 1) {
	$left{'chr0'} =~ m/(.*)/i; # this puts the base chromosome name ('chr.') into $1
	my $tmp = $1;
	$tmp =~ tr/A-Z/a-z/;
	$tmp =~ s/rc_//i;
	my $match = checkMatch ($left{'chr0'}, $right{'chr0'}, $left{'coord'}, $right{'coord'}, $reference{$tmp}, $offset, $distance);

	if($match) {
	    if ($left{'chr0'} =~ /^RC_/i) {

                $l_seqname   = $tmp ;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $reference{$tmp} + 1 - $left{'coord'};
                $l_end       = $reference{$tmp} + 1 - $left{'coord'} + $readsize;
                $l_score     = 1;
                $l_strand    = q{-};
                $l_attribute = 'target=' . substr ($reference{"$tmp-rc"}, $left{'coord'} - 3, $readsize + 4);

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $right{'coord'};
                $r_end       = $right{'coord'} + $readsize;
                $r_score     = 1;
                $r_strand    = q{+};
                $r_attribute = 'target=' . substr ($reference{"$tmp-seq"}, $right{'coord'} - 3, $readsize + 4);

                ### TO DELETE
                print join("\t",
			   $tmp,
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   $reference{$tmp} + 1 - $left{'coord'},
			   $reference{$tmp} + 1 - $left{'coord'} + $readsize,
			   1,
			   "-",
			   ".",
			   "target=" . substr($reference{"$tmp-rc"}, $left{'coord'}-3, $readsize+4)) . "\n";
		print join("\t",
			   $tmp,
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   $right{'coord'},
			   $right{'coord'}+$readsize,
			   1,
			   "+",
			   ".",
			   "target=" . substr($reference{"$tmp-seq"}, $right{'coord'}-3, $readsize+4)) . "\n";
                ### TO DELETE
	    }
	    else {

                $l_seqname   = $tmp ;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $left{'coord'};
                $l_end       = $left{'coord'} + $readsize;
                $l_score     = 1;
                $l_strand    = q{+};
                $l_attribute = 'target=' . substr ($reference{"$tmp-seq"}, $left{'coord'} - 3, $readsize + 4);

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $reference{$tmp} + 1 - $right{'coord'};
                $r_end       = $reference{$tmp} + 1 - $right{'coord'} + $readsize;
                $r_score     = 1;
                $r_strand    = q{-};
                $r_attribute = 'target=' . substr ($reference{"$tmp-rc"}, $right{'coord'} - 3, $readsize + 4);

                ### TO DELETE
                print join("\t",
			   $tmp,
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   $left{'coord'},
			   $left{'coord'} + $readsize,
			   1,
			   "+",
			   ".",
			   "target=" . substr($reference{"$tmp-seq"}, $left{'coord'}-3, $readsize+4)) . "\n";
		print join("\t",
			   $tmp,
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   $reference{$tmp} + 1 - $right{'coord'},
			   $reference{$tmp} + 1 - $right{'coord'} + $readsize,
			   1,
			   "-",
			   ".",
			   "target=" . substr($reference{"$tmp-rc"}, $right{'coord'}-3, $readsize+4)) . "\n";
                ### TO DELETE
	    }
	}
	else { # if not match
	    if ($nomatches) {

                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $r_feature   = $right{'line'} . q{:} . $rsequence;

                ### TO DELETE
		print join("\t",
			   ".",
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   ".",
			   ".",
			   0,
			   ".",
			   ".",
			   ".") . "\n";
		print join("\t",
			   ".",
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   ".",
			   ".",
			   0,
			   ".",
			   ".",
			   ".") . "\n";
                ### TO DELETE
	    }
	}
    }
    ##### Unique matches on both ends #####



    ##### NM - U0 #####
    ##### No match on one end, 1 match on other #####
    if ($lmatch == 0 && $rmatch == 1) {
        $right{'chr0'}=~m/(.*)/i; # this puts the base chromosome name ('chr.') into $1
        my $tmp=$1;
        $tmp=~tr/A-Z/a-z/;
	$tmp=~s/rc_//i;

	if ($nomatches) {

            $l_feature   = $left{'line'} . q{:} . $lsequence;

            ### TO DELETE
	    print join("\t",
		       ".",
		       "pair_align",
		       $left{'line'} . ":" . $lsequence,
		       ".",
		       ".",
		       0,
		       ".",
		       ".",
		       ".") . "\n";
            ### TO DELETE
	}

        if ($right{'chr0'} =~ m/^RC_/i) { # checks that this maps to reverse strand

            $r_seqname   = $tmp;
            $r_feature   = $right{'line'} . q{:} . $rsequence;
            $r_start     = $reference{$tmp} + 1 - $right{'coord'};
            $r_end       = $reference{$tmp} + 1 - $right{'coord'} + $readsize;
            $r_score     = 1;
            $r_strand    = q{-};
            $r_attribute = 'target=' . substr ($reference{"$tmp-rc"}, $right{'coord'} - 3, $readsize + 4);

            ### TO DELETE
	    print join("\t",
		       $tmp,
		       "pair_align",
		       $right{'line'} . ":" . $rsequence,
		       $reference{$tmp} + 1 - $right{'coord'},
		       $reference{$tmp} + 1 - $right{'coord'} + $readsize,
		       1,
		       "-",
		       ".",
		       "target=" . substr($reference{"$tmp-rc"}, $right{'coord'}-3, $readsize+4)) . "\n";
            ### TO DELETE
        }
        else {

            $r_seqname   = $tmp;
            $r_feature   = $right{'line'} . q{:} . $rsequence;
            $r_start     = $right{'coord'};
            $r_end       = $right{'coord'} + $readsize;
            $r_score     = 1;
            $r_strand    = q{+};
            $r_attribute = 'target=' . substr ($reference{"$tmp-seq"}, $right{'coord'} - 3, $readsize + 4);

            ### TO DELETE
	    print join("\t",
		       $tmp,
		       "pair_align",
		       $right{'line'} . ":" . $rsequence,
		       $right{'coord'},
		       $right{'coord'} + $readsize,
		       1,
		       "+",
		       ".",
		       "target=" . substr($reference{"$tmp-seq"}, $right{'coord'}-3, $readsize+4)) . "\n";
            ### TO DELETE
        }
    }

    if ($lmatch == 1 && $rmatch == 0) {
        $left{'chr0'} =~ m/(.*)/i; # this puts the base chromosome name ('chr.') into $1
        my $tmp = $1;
        $tmp =~ tr/A-Z/a-z/;
	$tmp =~ s/rc_//i;

	if($left{'chr0'}=~m/^RC_/i) {

            $l_seqname   = $tmp;
            $l_feature   = $left{'line'} . q{:} . $lsequence;
            $l_start     = $reference{$tmp} + 1 - $left{'coord'};
            $l_end       = $reference{$tmp} + 1 - $left{'coord'} + $readsize;
            $l_score     = 1;
            $l_strand    = q{-};
            $l_attribute = 'target=' . substr ($reference{"$tmp-rc"}, $left{'coord'} - 3, $readsize + 4);

            ### TO DELETE
            print join("\t",
		       $tmp,
		       "pair_align",
		       $left{'line'} . ":" . $lsequence,
		       $reference{$tmp} + 1 - $left{'coord'},
		       $reference{$tmp} + 1 - $left{'coord'} + $readsize,
		       1,
		       "-",
		       ".",
		       "target=" . substr($reference{"$tmp-rc"}, $left{'coord'}-3, $readsize+4)) . "\n";
            ### TO DELETE
        }
        else {

            $l_seqname   = $tmp;
            $l_feature   = $left{'line'} . q{:} . $lsequence;
            $l_start     = $left{'coord'};
            $l_end       = $left{'coord'} + $readsize;
            $l_score     = 1;
            $l_strand    = q{+};
            $l_attribute = 'target=' . substr ($reference{"$tmp-seq"}, $left{'coord'} - 3, $readsize + 4);

            ### TO DELETE
            print join("\t",
		       $tmp,
		       "pair_align",
		       $left{'line'} . ":" . $lsequence,
		       $left{'coord'},
		       $left{'coord'} + $readsize,
		       1,
		       "+",
		       ".",
		       "target=" . substr($reference{"$tmp-seq"}, $left{'coord'}-3, $readsize+4)) . "\n";
            ### TO DELETE
        }

	if ($nomatches) {

            $r_feature   = $right{'line'} . q{:} . $rsequence;

            ### TO DELETE
	    print join("\t",
		       ".",
		       "pair_align",
		       $right{'line'} . ":" . $rsequence,
		       ".",
		       ".",
		       0,
		       ".",
		       ".",
		       ".") . "\n";
            ### TO DELETE
	}
    }
    ##### No match on one end, 1 match on other #####



    ##### NM - R0 #####
    ##### No match on one end, multiple matches on other #####
    if ($repeats) {
	if ($lmatch == 0 && $rmatch > 1) {
	    if ($nomatches) {

                $l_feature   = $left{'line'} . q{:} . $lsequence;

                ### TO DELETE
                print join("\t",
			   ".",
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   ".",
			   ".",
			   0,
			   ".",
			   ".",
			   ".") . "\n";
                ### TO DELETE
	    }

            $r_feature   = $right{'line'} . q{:} . $rsequence;
            $r_score     = $rmatch;
            $r_attribute = 'targets=';

            ### TO DELETE
            print join("\t",
		       ".",
		       "pair_align",
		       $right{'line'} . ":" . $rsequence,
		       ".",
		       ".",
		       $rmatch,
		       ".",
		       ".",
		       "targets=");
            ### TO DELETE

	    for my $i (0 .. $rmatch - 1) {
		$right{"chr$i"} =~ m/(.*)/i; # this puts the base chromosome name ('chr.') into $1
		my $tmp = $1;
		$tmp =~ tr/A-Z/a-z/;
		$tmp =~ s/rc_//i;

		if ($right{"chr$i"} =~ m/^RC_/i) {

                    $r_attribute .= join(q{:},
                                         $right{"chr$i"},
                                         $reference{$tmp} + 1 - $right{"coord$i"},
                                         substr ($reference{"$tmp-rc"}, $right{"coord$i"} - 3, $readsize + 4)
                                     );
                    ### TO DELETE
		    print join(":",
			       $right{"chr$i"},
			       $reference{$tmp} + 1 - $right{"coord$i"},
			       substr($reference{"$tmp-rc"}, $right{"coord$i"}-3, $readsize+4));
                    ### TO DELETE
		}
		else {

                    $r_attribute .= join(q{:},
                                         $right{"chr$i"},
                                         $right{"coord$i"},
                                         substr ($reference{"$tmp-seq"}, $right{"coord$i"} - 3, $readsize + 4)
                                     );

                    ### TO DELETE
                    print join(":",
			       $right{"chr$i"},
			       $right{"coord$i"},
			       substr($reference{"$tmp-seq"}, $right{"coord$i"}-3, $readsize+4));
                    ### TO DELETE
		}

                $r_attribute .= q{,};

                ### TO DELETE
                print ",";
                ### TO DELETE
	    }

            $r_attribute .= "\b\n";

            ### TO DELETE
            print "\b\n";
            ### TO DELETE
	}

	if ($lmatch > 1 && $rmatch == 0) {

            $l_feature   = $left{'line'} . q{:} . $lsequence;
            $l_score     = $lmatch;
            $l_attribute = 'targets=';

            ### TO DELETE
	    print join("\t",
		       ".",
		       "pair_align",
		       $left{'line'} . ":" . $lsequence,
		       ".",
		       ".",
		       $lmatch,
		       ".",
		       ".",
		       "targets=");
            ### TO DELETE

	    for my $i (0..$lmatch - 1) {
		$left{"chr$i"}=~m/(.*)/i; # this puts the base chromosome name ('chr.') into $1
		my $tmp=$1;
		$tmp=~tr/A-Z/a-z/;
		$tmp=~s/rc_//i;

		if($left{"chr$i"}=~m/^RC_/i) {

                    $l_attribute .= join(q{:},
                                         $left{"chr$i"},
                                         $reference{$tmp} + 1 - $left{"coord$i"},
                                         substr ($reference{"$tmp-rc"}, $left{"coord$i"} - 3, $readsize + 4)
                                     );

                    ### TO DELETE
                    print join(":",
			       $left{"chr$i"},
			       $reference{$tmp} + 1 - $left{"coord$i"},
			       substr($reference{"$tmp-rc"}, $left{"coord$i"}-3, $readsize+4));
                    ### TO DELETE
		}
		else {

                    $l_attribute .= join(q{:},
                                         $left{"chr$i"},
                                         $left{"coord$i"},
                                         substr ($reference{"$tmp-seq"}, $left{"coord$i"} - 3, $readsize + 4)
                                     );

                    ### TO DELETE
                    print join(":",
			       $left{"chr$i"},
			       $left{"coord$i"},
			       substr($reference{"$tmp-seq"}, $left{"coord$i"}-3, $readsize+4));
                    ### TO DELETE
		}

                $l_attribute .= q{,};

                ### TO DELETE
		print ",";
                ### TO DELETE
	    }

            $l_attribute .= "\b\n";

            ### TO DELETE
	    print "\b\n";
            ### TO DELETE

	    if ($nomatches) {

                $r_feature   = $right{'line'} . q{:} . $rsequence;

                ### TO DELETE
                print join("\t",
			   ".",
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   ".",
			   ".",
			   0,
			   ".",
			   ".",
			   ".") . "\n";
                ### TO DELETE
	    }
	}
    }
    ##### No match on one end, multiple matches on other #####



    ##### One match on one end, multiple matches on other end #####
    if ($lmatch == 1 && $rmatch > 1) {
	$left{'chr0'} =~ m/(.*)/i; # this puts the base chromosome name ('chr.') into $1
	my $tmp = $1;
	$tmp =~ tr/A-Z/a-z/;
	$tmp =~ s/rc_//i;

 	my ($bestcoord, $beststrand);
	my $score = 100000000;
 	for my $i (0..$rmatch - 1) {
	    my $tmpscore = checkMultMatch ($left{'chr0'}, $right{"chr$i"}, $left{'coord'}, $right{"coord$i"}, $reference{$tmp}, $offset, $distance);
	    if($tmpscore == -1) {next}
	    if($tmpscore < $score) {
		$score = $tmpscore;
		$bestcoord = $right{"coord$i"};
		$beststrand = $right{"chr$i"};
	    }
	    else {next}
 	}

	$beststrand =~ m/(.*)/i; # this puts the base chromosome name ('chr.') into $1
	$tmp = $1;
	$tmp =~ tr/A-Z/a-z/;
	$tmp =~ s/rc_//i;
	if ($beststrand =~ m/^RC_/i) {
	    if ($score <= $offset + $distance) {

                $l_seqname   = $tmp;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $left{'coord'};
                $l_end       = $left{'coord'} + $readsize;
                $l_score     = 1;
                $l_strand    = q{+};
                $l_attribute = 'target=' . substr ($reference{"$tmp-seq"}, $left{'coord'} - 3, $readsize + 4);

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $reference{$tmp} + 1 - $bestcoord;
                $r_end       = $reference{$tmp} + 1 - $bestcoord + $readsize;
                $r_score     = 1;
                $r_strand    = q{-};
                $r_attribute = 'target=' . substr ($reference{"$tmp-rc"}, $bestcoord - 3, $readsize + 4);

                ### TO DELETE
                print join("\t",
			   $tmp,
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   $left{'coord'},
			   $left{'coord'} + $readsize,
			   1,
			   "+",
			   ".",
			   "target=" . substr($reference{"$tmp-seq"}, $left{'coord'}-3, $readsize+4)) . "\n";
		print join("\t",
			   $tmp,
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   $reference{$tmp} + 1 - $bestcoord,
			   $reference{$tmp} + 1 - $bestcoord + $readsize,
			   1,
			   "-",
			   ".",
			   "target=" . substr($reference{"$tmp-rc"}, $bestcoord-3, $readsize+4)) . "\n";
                ### TO DELETE
	    }
	    else {

                $l_seqname   = $tmp;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $reference{$tmp} + 1 - $left{'coord'};
                $l_end       = $reference{$tmp} + 1 - $left{'coord'} + $readsize;
                $l_score     = 1;
                $l_strand    = q{-};
                $l_attribute = 'target=' . substr ($reference{"$tmp-rc"}, $left{'coord'} - 3, $readsize + 4);

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $bestcoord;
                $r_end       = $bestcoord + $readsize;
                $r_score     = 1;
                $r_strand    = q{+};
                $r_attribute = 'target=' . substr ($reference{"$tmp-seq"}, $bestcoord - 3, $readsize + 4);

                ### TO DELETE
                print join("\t",
			   $tmp,
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   $reference{$tmp} + 1 - $left{'coord'},
			   $reference{$tmp} + 1 - $left{'coord'} + $readsize,
			   1,
			   "-",
			   ".",
			   "target=" . substr($reference{"$tmp-rc"}, $left{'coord'}-3, $readsize+4)) . "\n";
		print join("\t",
			   $tmp,
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   $bestcoord,
			   $bestcoord + $readsize,
			   1,
			   "+",
			   ".",
			   "target=" . substr($reference{"$tmp-seq"}, $bestcoord-3, $readsize+4)) . "\n";
                ### TO DELETE
	    }
	}
	else {
	    if ($nomatches) {

                $l_feature  = $left{'line'} . q{:} . $lsequence;
                $r_feature  = $right{'line'} . q{:} . $rsequence;

                ### TO DELETE
		print join("\t",
			   ".",
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   ".",
			   ".",
			   0,
			   ".",
			   ".",
			   ".") . "\n";

		print join("\t",
			   ".",
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   ".",
			   ".",
			   0,
			   ".",
			   ".",
			   ".") . "\n";
                ### TO DELETE
	    }
	}
    }

    if ($lmatch > 1 && $rmatch == 1) {
	$right{'chr0'} =~ m/(.*)/i; # this puts the base chromosome name ('chr.') into $1
	my $tmp = $1;
	$tmp =~ tr/A-Z/a-z/;
	$tmp =~ s/rc_//i;

 	my ($bestcoord, $beststrand);
	my $score = 100000000;
 	for my $i (0..$lmatch - 1) {
	    my $tmpscore = checkMultMatch ($left{"chr$i"}, $right{'chr0'}, $left{"coord$i"}, $right{'coord'}, $reference{$tmp}, $offset, $distance);
	    if ($tmpscore == -1) {next}
	    if ($tmpscore < $score) {
		$score = $tmpscore;
		$bestcoord = $left{"coord$i"};
		$beststrand = $left{"chr$i"};
	    }
            else {next}
 	}

	$beststrand =~ m/(.*)/i; # this puts the base chromosome name ('chr.') into $1
	$tmp = $1;
	$tmp =~ tr/A-Z/a-z/;
	$tmp =~ s/rc_//i;

	if ($beststrand =~ m/^RC_/i) {
	    if ($score <= $offset + $distance) {

                $l_seqname   = $tmp;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $reference{$tmp} + 1 - $bestcoord;
                $l_end       = $reference{$tmp} + 1 - $bestcoord + $readsize;
                $l_score     = 1;
                $l_strand    = q{-};
                $l_attribute = 'target=' . substr ($reference{"$tmp-rc"}, $bestcoord - 3, $readsize + 4);

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $right{'coord'};
                $r_end       = $right{'coord'} + $readsize;
                $r_score     = 1;
                $r_strand    = q{+};
                $r_attribute = 'target=' . substr ($reference{"$tmp-seq"}, $right{'coord'} - 3, $readsize + 4);

                ### TO DELETE
		print join("\t",
			   $tmp,
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   $reference{$tmp} + 1 - $bestcoord,
			   $reference{$tmp} + 1 - $bestcoord + $readsize,
			   1,
			   "-",
			   ".",
			   "target=" . substr($reference{"$tmp-rc"}, $bestcoord-3, $readsize+4)) . "\n";
		print join("\t",
			   $tmp,
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   $right{'coord'},
			   $right{'coord'} + $readsize,
			   1,
			   "+",
			   ".",
			   "target=" . substr($reference{"$tmp-seq"}, $right{'coord'}-3, $readsize+4)) . "\n";
                ### TO DELETE
	    }
	    else {

                $l_seqname   = $tmp;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $bestcoord;
                $l_end       = $bestcoord + $readsize;
                $l_score     = 1;
                $l_strand    = q{+};
                $l_attribute = 'target=' . substr ($reference{"$tmp-seq"}, $bestcoord - 3, $readsize + 4);

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $reference{$tmp} + 1 - $right{'coord'};
                $r_end       = $reference{$tmp} + 1 - $right{'coord'} + $readsize;
                $r_score     = 1;
                $r_strand    = q{-};
                $r_attribute = 'target=' . substr ($reference{"$tmp-rc"}, $right{'coord'} - 3, $readsize + 4);

                ### TO DELETE
                print join("\t",
			   $tmp,
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   $bestcoord,
			   $bestcoord + $readsize,
			   1,
			   "+",
			   ".",
			   "target=" . substr($reference{"$tmp-seq"}, $bestcoord-3, $readsize+4)) . "\n";
		print join("\t",
			   $tmp,
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   $reference{$tmp} + 1 - $right{'coord'},
			   $reference{$tmp} + 1 - $right{'coord'} + $readsize,
			   1,
			   "-",
			   ".",
			   "target=" . substr($reference{"$tmp-rc"}, $right{'coord'}-3, $readsize+4)) . "\n";
                ### TO DELETE
	    }
	}
	else {
	    if ($nomatches) {

                $l_feature = $left{'line'} . ":" . $lsequence;
                $r_feature = $right{'line'} . ":" . $rsequence;

                ### TO DELETE
		print join("\t",
			   ".",
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   ".",
			   ".",
			   0,
			   ".",
			   ".",
			   ".") . "\n";

		print join("\t",
			   ".",
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   ".",
			   ".",
			   0,
			   ".",
			   ".",
			   ".") . "\n";
                ### TO DELETE
	    }
	}
    }
    ##### One match on one end, multiple matches on other end #####



    ##### R0 - R0 #####
    ##### Multiple matches on both ends #####
    if ($lmatch > 1 && $rmatch > 1) {
	my ($lbestcoord, $rbestcoord, $lbeststrand, $rbeststrand);
	my $bestscore = 100000000;

	# loops through each possible target on the left
	for my $i (0..$lmatch - 1) {

	    $left{"chr$i"} =~ m/(.*)/i; # this puts the base chromosome name ('chr.') into $1
	    my $tmp = $1;
	    $tmp =~ tr/A-Z/a-z/;
	    $tmp =~ s/rc_//i;

	    my ($bestcoord, $beststrand);
	    my $score = 100000000;

	    # loops through each possible target on the right
	    for my $j (0..$rmatch - 1) {

		my $tmpscore = checkMultMatch($left{"chr$i"}, $right{"chr$j"}, $left{"coord$i"}, $right{"coord$j"}, $reference{$tmp}, $offset, $distance);

		if($tmpscore == -1) {next}
		if($tmpscore < $score) {
		    $score = $tmpscore;
		    $bestcoord = $right{"coord$j"};
		    $beststrand = $right{"chr$j"};
		}
	    }

	    # if no good matches on right side, skip this left sequence
	    if($bestcoord == -1 && $beststrand == -1) {next}

	    if($score < $bestscore) {
		$bestscore = $score;
		$lbestcoord = $left{"coord$i"};
		$lbeststrand = $left{"chr$i"};
		$rbestcoord = $bestcoord;
		$rbeststrand = $beststrand;
	    }
	}

	$lbeststrand =~ m/(.*)/i; # this puts the base chromosome name ('chr.') into $1
	my $tmp = $1;
	$tmp =~ tr/A-Z/a-z/;
	$tmp =~ s/rc_//i;

	if ($lbeststrand =~ m/^RC_/i) {
	    if ($bestscore <= $offset + $distance) {

                $l_seqname   = $tmp;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $reference{$tmp} + 1 - $lbestcoord;
                $l_end       = $reference{$tmp} + 1 - $lbestcoord + $readsize;
                $l_score     = 1;
                $l_strand    = q{-};
                $l_attribute = 'target=' . substr ($reference{"$tmp-rc"}, $lbestcoord - 3, $readsize + 4);

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $rbestcoord;
                $r_end       = $rbestcoord + $readsize;
                $r_score     = 1;
                $r_strand    = q{+};
                $r_attribute = 'target=' . substr ($reference{"$tmp-seq"}, $rbestcoord - 3, $readsize + 4);

                ### TO DELETE
                print join("\t",
			   $tmp,
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   $reference{$tmp} + 1 - $lbestcoord,
			   $reference{$tmp} + 1 - $lbestcoord + $readsize,
			   1,
			   "-",
			   ".",
			   "target=" . substr($reference{"$tmp-rc"}, $lbestcoord-3, $readsize+4)) . "\n";
		print join("\t",
			   $tmp,
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   $rbestcoord,
			   $rbestcoord + $readsize,
			   1,
			   "+",
			   ".",
			   "target=" . substr($reference{"$tmp-seq"}, $rbestcoord-3, $readsize+4)) . "\n";
                ### TO DELETE
	    }
	    else {

                $l_seqname   = $tmp;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $lbestcoord;
                $l_end       = $lbestcoord + $readsize;
                $l_score     = 1;
                $l_strand    = q{+};
                $l_attribute = 'target=' . substr ($reference{"$tmp-seq"}, $lbestcoord - 3, $readsize + 4);

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $reference{$tmp} + 1 - $rbestcoord;
                $r_end       = $reference{$tmp} + 1 - $rbestcoord + $readsize;
                $r_score     = 1;
                $r_strand    = q{-};
                $r_attribute = 'target=' . substr ($reference{"$tmp-rc"}, $rbestcoord - 3, $readsize + 4);

                ### TO DELETE
		print join("\t",
			   $tmp,
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   $lbestcoord,
			   $lbestcoord + $readsize,
			   1,
			   "+",
			   ".",
			   "target=" . substr($reference{"$tmp-seq"}, $lbestcoord-3, $readsize+4)) . "\n";
		print join("\t",
			   $tmp,
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   $reference{$tmp} + 1 - $rbestcoord,
			   $reference{$tmp} + 1 - $rbestcoord + $readsize,
			   1,
			   "-",
			   ".",
			   "target=" . substr($reference{"$tmp-rc"}, $rbestcoord-3, $readsize+4)) . "\n";
                ### TO DELETE
	    }
	}
	else {
	    if ($nomatches) {

                $l_feature = $left{'line'} . ":" . $lsequence;
                $r_feature = $right{'line'} . ":" . $rsequence;

                ### TO DELETE
                print join("\t",
			   ".",
			   "pair_align",
			   $left{'line'} . ":" . $lsequence,
			   ".",
			   ".",
			   0,
			   ".",
			   ".",
			   ".") . "\n";
		print join("\t",
			   ".",
			   "pair_align",
			   $right{'line'} . ":" . $rsequence,
			   ".",
			   ".",
			   0,
			   ".",
			   ".",
			   ".") . "\n";
                ### TO DELETE
	    }
	}
    }
    ##### Multiple matches on both ends #####
}
# end for loop through all sequences
print STDERR "done!\n";
close STDOUT;
close $LEFT;
close $RIGHT;
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
    if($lchr=~m/^([^R][^C][^_])/i) {
	if(!$rchr=~m/^RC_$1/i) {$match=0;}
    }
    else {
	if($rchr=~m/^RC_(.*)/i) {
	    if($lchr=~m/^$1/i) {$match=0;}
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
    if($lchr=~m/^([^R][^C][^_])/i) {
	if(!$rchr=~m/^RC_$1/i) {$score=-1;}
    }
    else {
	if($rchr=~m/^RC_(.*)/i) {
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

    if($hash{'matches'}>1) {
	my @all=split(/,/, $line[3]);
	$hash{'matches'}=scalar(@all);
    }

    if($hash{'matches'}==1) {
	my $tmp=(split ':', $line[3])[0];
	$tmp=~s/\n//g;
	$hash{'chr0'}=$tmp;
	$tmp=(split ':', $line[3])[1];
	$tmp=~s/[A-Z][0-9]$//i;
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


sub usage {
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
\t[--nomatches]\tPrint reads with no matches
\t[--quiet]\tDon't output perl warnings
\t[--verbose]\tOutput detailed perl diagnostic messages
\t[--usage]\tPrints this
Takes a paired ends sequencing output (needs to be preprocessed to fit above) 
and compares each pair to maximize unique alignment matches to a genome.\n";
	exit 1;
    }
}
