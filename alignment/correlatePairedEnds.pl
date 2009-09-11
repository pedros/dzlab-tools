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
use Data::Dumper;
#use Smart::Comments '###';
#use Text::LevenshteinXS qw(distance);
use List::Util qw(min);

# Globals, passed as command line options
my $leftendfile;           # left end sequences file in eland 3 format
my $rightendfile;          # right end sequences file in eland 3 format
my $referencefile;         # unmodified reference genome file in fasta format
my $skip_nm        = 0;
my $distance       = 300;  # maximum (random) range distance beyond offset.
my $offset         = 0;    # minimum distance between left and right sequences of a pair.
my $output         = q{-}; # output mode. standard output by default
my $readsize       = 45;   # length of each read
my $usage          = 0;    # print usage and exit
my $trust_dash_2   = 0;
my $single_ends    = 0;

# Initial check of command line parameters
usage();

# Parse command line arguments
my $result = GetOptions(
    "left|l=s"         => \$leftendfile,
    "right|r=s"        => \$rightendfile,
    "reference|f=s"    => \$referencefile,
    'skip-nm|k'        => \$skip_nm,
    "distance|d=i"     => \$distance,
    "offset|t=i"       => \$offset,
    'single-ends|1=i'  => \$single_ends,
    'trust-dash-2|2=i' => \$trust_dash_2,
    "output|o:s"       => \$output,
    "readsize|s:i"     => \$readsize,
    "verbose|v"        => sub { enable diagnostics },
    "quiet|q"          => sub { no warnings  },
    "help|h"           => \&usage
);

# holds name of chromosomes as keys and length of chromosomes in bp as values
my %reference = ();

# begins limited scope for extracting names and lengths of chromosomes in reference file
{
    # reads in the reference genome file into @fastaseq
    open my $REF, '<', "$referencefile" or croak "Can't open file: $referencefile";
    my @fastaseq = <$REF>;
    close $REF;

    # find and store indices for each chromosome change and corresponding descriptions
    my ( @idx, @dsc ) = ();
    for my $i ( 0 .. @fastaseq - 1 ) { ### Indexing $referencefile...  % done
        if ( $fastaseq[$i] =~ m/^>/ ) {
            $fastaseq[$i] =~ s/>//g;
            $fastaseq[$i] = ( split /\s/, "$fastaseq[$i]" )[0];
            $fastaseq[$i] =~ tr/A-Z/a-z/;
            push @idx, $i;
            push @dsc, $fastaseq[$i];
        }
    }

    # gets and saves each chromosome's sequence and reverse complemented sequence
    for my $j ( 0 .. @idx - 1 ) { ### Loading $referencefile into memory...  % done
        my $line;
        if ( $j == scalar @idx - 1 ) {
            $line = join( q{}, 'NN', @fastaseq[ $idx[$j] + 1 .. @fastaseq - 1] , 'NN');
        }
        else {
            $line = join( q{}, 'NN', @fastaseq[ $idx[$j] + 1 .. $idx[$j + 1] - 1] , 'NN');
        }
        $line =~ s/[\n\r]//g;
        $reference{ $dsc[$j] }     = (length $line) - 4;
        $reference{"$dsc[$j]-seq"} = $line;
        $reference{"$dsc[$j]-rc"}  = reverseComp ($line);
    }
}

# redirects STDOUT to file if specified by user
if ( $output ne q{-} ) {
    open STDOUT, '>', "$output";
}

# opens sequence files
open my $LEFT,  '<', $leftendfile
or croak "Can't open file: $leftendfile";
open my $RIGHT, '<', $rightendfile
or croak "Can't open file: $rightendfile";

while (my $leftend = <$LEFT>) {
    # reads single sequence from each file
    my $rightend = <$RIGHT>;

    $leftend =~ s/[\r\n]//g;
    $rightend =~ s/[\r\n]//g;

    # parses each line into hash
    my %left  = %{ parseEland3Line($leftend) };
    my %right = %{ parseEland3Line($rightend) };

    # gets number of matches and sequence from hash
    my $lmatch    = $left{'matches'};
    my $rmatch    = $right{'matches'};
    my $lsequence = $left{'sequence'};
    my $rsequence = $right{'sequence'};
    my $lreadsize = length $lsequence;
    my $rreadsize = length $rsequence;

    # initializes each gff field to default values
    my (
        $l_seqname, $l_source, $l_feature, $l_start, $l_end,
        $l_score,   $l_strand, $l_frame,   $l_attribute
    ) = ( q{.}, 'pair_align', q{.}, 0, 0, 0, q{.}, q{.}, q{.} );
    my (
        $r_seqname, $r_source, $r_feature, $r_start, $r_end,
        $r_score,   $r_strand, $r_frame,   $r_attribute
    ) = ( q{.}, 'pair_align', q{.}, 0, 0, 0, q{.}, q{.}, q{.} );



    ##### START POSSIBLE CASES HERE #####

    ##### No matches on either end #####
    if ( $lmatch == 0 && $rmatch == 0 ) {
        $l_feature = $left{'line'} . q{:} . $lsequence;
        $r_feature = $right{'line'} . q{:} . $rsequence;
        $l_source  = 'NM/NM';
        $r_source  = 'NM/NM';
    }
    ##### No matches on either end #####



    ##### Unique matches on both ends #####
    if ( $lmatch == 1 && $rmatch == 1 ) {

        $l_source  = 'U/U';
        $r_source  = 'U/U';

        # gets chromosome name into $tmp
        $left{'chr0'} =~ m/(.*)/i;
        my $tmp = $1;
        $tmp =~ tr/A-Z/a-z/;
        $tmp =~ s/rc_//i;

        # checks that both sequences match
        my $match =
            check_match( $left{'chr0'}, $right{'chr0'}, $left{'coord'},
                        $right{'coord'}, $reference{$tmp}, $offset, $distance, $single_ends, $lreadsize, $rreadsize);

        # if both sequences match
        if ($match) {
            if ( $left{'chr0'} =~ /^RC_/i ) { # if left sequence maps to reverse strand
                $l_seqname = $tmp;
                $l_feature = $left{'line'} . q{:} . $lsequence;
                $l_start   = $left{'coord'};
                $l_end     = $left{'coord'} + $lreadsize - 1;
                $l_score   = 1;
                $l_strand  = q{-};
                $l_frame   = $left{mm};
                $l_attribute = 'target='
                    . substr(
                        $reference{"$tmp-rc"},
                        $left{'coord'} - 1,
                        $lreadsize + 4
                    );

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $right{'coord'};
                $r_end       = $right{'coord'} + $rreadsize - 1;
                $r_score     = 1;
                $r_frame   = $right{mm};

                if ($single_ends) {
                    $r_strand    = q{-};
                    $r_attribute = 'target='
                    . substr(
                        $reference{"$tmp-rc"},
                        $right{'coord'} - 1,
                        $rreadsize + 4
                    );
                }
                else {
                    $r_strand    = q{+};
                    $r_attribute = 'target='
                    . substr(
                        $reference{"$tmp-seq"},
                        $right{'coord'} - 1,
                        $rreadsize + 4
                    );
                }
            }
            else {         # if left sequence maps to forward strand
                $l_seqname   = $tmp;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $left{'coord'};
                $l_end       = $left{'coord'} + $lreadsize - 1;
                $l_score     = 1;
                $l_strand    = q{+};
                $l_frame     = $left{mm};
                $l_attribute = 'target='
                    . substr(
                        $reference{"$tmp-seq"},
                        $left{'coord'} - 1,
                        $lreadsize + 4
                    );

                $r_seqname = $tmp;
                $r_feature = $right{'line'} . q{:} . $rsequence;
                $r_start   = $right{'coord'};
                $r_end     = $right{'coord'} + $rreadsize - 1;
                $r_score   = 1;
                $r_frame     = $right{mm};

                if ($single_ends) {
                    $r_strand  = q{+};
                    $r_attribute = 'target='
                    . substr(
                        $reference{"$tmp-seq"},
                        $right{'coord'} - 1,
                        $rreadsize + 4
                    );
                }
                else {
                    $r_strand  = q{-};
                    $r_attribute = 'target='
                    . substr(
                        $reference{"$tmp-rc"},
                        $right{'coord'} - 1,
                        $rreadsize + 4
                    );
                }
            }
        } # no match, trust /1 read
        else {
            if ( $left{'chr0'} =~ /^RC_/i ) { # if left sequence maps to reverse strand
                $l_seqname = $tmp;
                $l_feature = $left{'line'} . q{:} . $lsequence;
                $l_start   = $left{'coord'};
                $l_end     = $left{'coord'} + $lreadsize - 1;
                $l_score   = 1;
                $l_strand  = q{-};
                $l_frame     = $left{mm};
                $l_attribute = 'target='
                    . substr(
                        $reference{"$tmp-rc"},
                        $left{'coord'} - 1,
                        $lreadsize + 4
                    );
            }
            else {         # if left sequence maps to forward strand
                $l_seqname   = $tmp;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $left{'coord'};
                $l_end       = $left{'coord'} + $lreadsize - 1;
                $l_score     = 1;
                $l_strand    = q{+};
                $l_frame     = $left{mm};
                $l_attribute = 'target='
                . substr(
                    $reference{"$tmp-seq"},
                    $left{'coord'} - 1,
                    $lreadsize + 4
                );
            }
            $r_feature = $right{'line'} . q{:} . $rsequence;
        }
    }

    ##### Unique matches on both ends #####



    ##### No match on one end, 1 match on other #####

    # 0 matches on left sequence
    if ( $lmatch == 0 && $rmatch == 1 ) {

        $l_source  = 'NM/U';
        $r_source  = 'NM/U';

        $l_feature = $left{'line'} . q{:} . $lsequence;
        $r_feature = $right{'line'} . q{:} . $rsequence;

        if ($trust_dash_2) {
            # gets chromosome name into $tmp
            $right{'chr0'} =~ m/(.*)/i;
            my $tmp = $1;
            $tmp =~ tr/A-Z/a-z/;
            $tmp =~ s/rc_//i;

            $l_feature = $left{'line'} . q{:} . $lsequence;

            if ( $right{'chr0'} =~ m/^RC_/i ) { # if right sequence maps to reverse strand
                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $right{'coord'};
                $r_end       = $right{'coord'} + $rreadsize - 1;
                $r_score     = 1;
                $r_strand    = q{-};
                $r_attribute = 'target='
                . substr(
                    $reference{"$tmp-rc"},
                    $right{'coord'} - 1,
                    $rreadsize + 4
                );
            }
            else {        # if right sequence maps to forward strand
                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $right{'coord'};
                $r_end       = $right{'coord'} + $rreadsize - 1;
                $r_score     = 1;
                $r_strand    = q{+};

                $r_attribute = 'target='
                . substr(
                    $reference{"$tmp-seq"},
                    $right{'coord'} - 1,
                    $rreadsize + 4
                );
            }
        }
    }

    # 0 matches on right sequence
    if ( $lmatch == 1 && $rmatch == 0 ) {

        $l_source  = 'U/NM';
        $r_source  = 'U/NM';

        # gets chromosome name into $tmp
        $left{'chr0'} =~ m/(.*)/i;
        my $tmp = $1;
        $tmp =~ tr/A-Z/a-z/;
        $tmp =~ s/rc_//i;

        if ( $left{'chr0'} =~ m/^RC_/i ) {  # if left sequence maps to reverse strand
            $l_seqname = $tmp;
            $l_feature = $left{'line'} . q{:} . $lsequence;
            $l_start   = $left{'coord'};
            $l_end     = $left{'coord'} + $lreadsize - 1;
            $l_score   = 1;
            $l_strand  = q{-};
            $l_frame   = $left{mm};
            $l_attribute =
            'target='
            . substr( $reference{"$tmp-rc"}, $left{'coord'} - 1,
                      $lreadsize + 4 );
        }
        else { # if left sequence maps to forward strand

            $l_seqname   = $tmp;
            $l_feature   = $left{'line'} . q{:} . $lsequence;
            $l_start     = $left{'coord'};
            $l_end       = $left{'coord'} + $lreadsize - 1;
            $l_score     = 1;
            $l_strand    = q{+};
            $l_frame   = $left{mm};
            $l_attribute = 'target='
            . substr(
                $reference{"$tmp-seq"},
                $left{'coord'} - 1,
                $lreadsize + 4
            );
        }
        $r_feature = $right{'line'} . q{:} . $rsequence;
    }
    ##### No match on one end, 1 match on other #####




    ##### No match on one end, multiple matches on other #####

    # multiple matches on right sequence
    if ( $lmatch == 0 && $rmatch > 1 ) {

        $l_source  = "NM/R";
        $r_source  = "NM/R";;

        $l_feature = $left{'line'} . q{:} . $lsequence;

        $r_feature   = $right{'line'} . q{:} . $rsequence;
        $r_score     = $rmatch;
        $r_attribute = 'targets=';

        # gets all possible matches and appends them to $r_attribute
        NM_MM_LOOP_RIGHT:
        for my $i ( 0 .. $rmatch - 1 ) {

            # gets the current chromosome
            $right{"chr$i"} =~ m/(.*)/i;
            my $tmp = $1;
            $tmp =~ tr/A-Z/a-z/;
            $tmp =~ s/rc_//i;

            if ( $right{"chr$i"} =~ m/^RC_/i ) { # if right sequence maps to reverse strand

                # appends chromosome name:coordinate:sequence to $r_attribute
                $r_attribute .= join(
                    q{:},
                    $right{"chr$i"},
                    $right{"coord$i"},
                    substr(
                        $reference{"$tmp-rc"},
                        $right{"coord$i"} - 1,
                        $rreadsize + 4
                    )
                );
            }
            else { # if right sequence maps to forward strand

                # appends chromosome name:coordinate:sequence to $r_attribute
                $r_attribute .= join(
                    q{:},
                    $right{"chr$i"},
                    $right{"coord$i"},
                    substr(
                        $reference{"$tmp-seq"},
                        $right{"coord$i"} - 1,
                        $rreadsize + 4
                    )
                );

            }
            $r_attribute .= q{,}; # appends sub-field separator comma
        }
        chop $r_attribute; # deletes very last comma and appends new line to $r_attribute
    }

    # multiple matches on left sequence
    if ( $lmatch > 1 && $rmatch == 0 ) {
        $l_source = "R/NM";
        $r_source = "R/NM";

        $l_feature   = $left{'line'} . q{:} . $lsequence;
        $l_score     = $lmatch;
        $l_attribute = 'targets=';

        # gets all possible matches and appends them to $r_attribute
        NM_MM_LOOP_LEFT:
        for my $i ( 0 .. $lmatch - 1 ) {
            $left{"chr$i"} =~ m/(.*)/i;
            my $tmp = $1;
            $tmp =~ tr/A-Z/a-z/;
            $tmp =~ s/rc_//i;

            if ( $left{"chr$i"} =~ m/^RC_/i ) { # if left sequence maps to reverse strand

                # appends chromosome name:coordinate:sequence to $l_attribute
                $l_attribute .= join(
                    q{:},
                    $left{"chr$i"},
                    $left{"coord$i"},
                    substr(
                        $reference{"$tmp-rc"},
                        $left{"coord$i"} - 1,
                        $lreadsize + 4
                    )
                );

            } else { # if left sequence maps to reverse strand

                # appends chromosome name:coordinate:sequence to $l_attribute
                $l_attribute .= join(
                    q{:},
                    $left{"chr$i"},
                    $left{"coord$i"},
                    substr(
                        $reference{"$tmp-seq"},
                        $left{"coord$i"} - 1,
                        $lreadsize + 4
                    )
                );
            }
            $l_attribute .= q{,}; # appends sub-field separator comma
        }
        chop $l_attribute; # deletes very last comma and appends new line to $r_attribute

        $r_feature = $right{'line'} . q{:} . $rsequence;
    }

    ##### No match on one end, multiple matches on other #####



    ##### One match on one end, multiple matches on other end #####

    # multiple matches on right sequence
    if ( $lmatch == 1 && $rmatch > 1 ) {

        $l_source  = 'U/R';
        $r_source  = 'U/R';

        # gets the left chromosome
        $left{'chr0'} =~ m/(.*)/i;
        my $tmp = $1;
        $tmp =~ tr/A-Z/a-z/;
        $tmp =~ s/rc_//i;

        if ( $left{'chr0'} =~ m/^RC_/i ) { # if left sequence maps to reverse strand
            $l_seqname   = $tmp;
            $l_feature   = $left{'line'} . q{:} . $lsequence;
            $l_start     = $left{'coord'};
            $l_end       = $left{'coord'} + $lreadsize - 1;
            $l_score     = 1;
            $l_strand    = q{-};
            $l_frame     = $left{mm};
            $l_attribute = 'target='
                . substr(
                    $reference{"$tmp-rc"},
                    $left{'coord'} - 1,
                    $lreadsize + 4
                );
        }
        else { # if left sequence maps to forward strand
            $l_seqname   = $tmp;
            $l_feature   = $left{'line'} . q{:} . $lsequence;
            $l_start     = $left{'coord'};
            $l_end       = $left{'coord'} + $lreadsize - 1;
            $l_score     = 1;
            $l_strand    = q{+};
            $l_frame     = $left{mm};
            $l_attribute = 'target='
                . substr(
                    $reference{"$tmp-seq"},
                    $left{'coord'} - 1,
                    $lreadsize + 4
                );
        }

        # loops through every possible match
        my ( $bestcoord, $beststrand, $bestmm ) = (-1, -1, -1);
        my $score = 100000000; # arbitrary large number; lower is better

        U0_MM_LOOP_RIGHT:
        for my $i ( 0 .. $rmatch - 1 ) {

            # checks best match
            my $tmpscore =
                check_mult_match( $left{'chr0'}, $right{"chr$i"}, $left{'coord'},
                                $right{"coord$i"}, $reference{$tmp}, $offset, $distance, $single_ends, $lreadsize, $rreadsize);

            # if $tmpscore is higher than current score of -1 discard potential match
            if ( $tmpscore == -1 or $tmpscore > $score) {next U0_MM_LOOP_RIGHT}

            # if $tmpscore flag is lower that current best score, keep
            if ( $tmpscore < $score ) {
                $score      = $tmpscore;
                $bestcoord  = $right{"coord$i"};
                $beststrand = $right{"chr$i"};
                $bestmm     = $right{"mm$i"};
            }
        }

        # gets the best right chromosome
        $beststrand =~ m/(.*)/i;
        $tmp = $1;
        $tmp =~ tr/A-Z/a-z/;
        $tmp =~ s/rc_//i;

        # if right sequence is distanced from left sequence within given range
        if ( $score >= $offset and $score <= $offset + $distance ) {

            # if right sequence maps to reverse strand
            if ( $beststrand =~ m/^RC_/i ) {

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $bestcoord;
                $r_end       = $bestcoord + $rreadsize - 1;
                $r_score     = 1;
                $r_strand    = q{-};
                $r_frame     = $bestmm,
                $r_attribute = 'target='
                    . substr(
                        $reference{"$tmp-rc"},
                        $bestcoord - 1,
                        $rreadsize + 4
                    );
            }
            else { # if right sequence maps to forward strand

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $bestcoord;
                $r_end       = $bestcoord + $rreadsize - 1;
                $r_score     = 1;
                $r_strand    = q{+};
                $r_frame     = $bestmm,
                $r_attribute = 'target='
                    . substr(
                        $reference{"$tmp-seq"},
                        $bestcoord - 1,
                        $rreadsize + 4
                    );

            }
        }
        else { # if right sequence is distanced from left sequence outside given range
            $r_feature = $right{'line'} . q{:} . $rsequence;
#             if ( $l_strand eq q{-} ) {

#                 my $right_target_sequence = substr $reference{"$l_seqname-rc"}, $left{'coord'} - 1 + $offset, $distance;
#                 $right_target_sequence =~ tr/Gg/Aa/;
#                 my @edit_distances = @{get_edit_distances ($rsequence, $right_target_sequence)};

#                 my $min_distance = min @edit_distances;
#                 my (@min_distance_indices) = grep { $edit_distances[$_] == $min_distance } 0..$#edit_distances;

#                 if (@min_distance_indices > 1) {
#                     $r_feature = $right{'line'} . q{:} . $rsequence;
#                 }
#                 else {
#                     $r_seqname   = $l_seqname;
#                     $r_feature   = $right{'line'} . q{:} . $rsequence;
#                     $r_start     = $left{'coord'} + $min_distance_indices[0];
#                     $r_end       = $left{'coord'} + $min_distance_indices[0] + $readsize;
#                     $r_score     = 1;
#                     $r_strand    = q{+};
#                     $r_attribute = 'target='
#                     . substr(
#                         $right_target_sequence,
#                         $min_distance_indices[0] - 2,
#                         $readsize + 4
#                     );
#                 }
#             }
#             else {

#                 my $right_target_sequence = substr $reference{"$l_seqname-seq"}, $left{'coord'} - 1 + $offset, $distance;
#                 $right_target_sequence =~ tr/Gg/Aa/;
#                 my @edit_distances = @{get_edit_distances ($rsequence, $right_target_sequence)};

#                 my $min_distance = min @edit_distances;
#                 my (@min_distance_indices) = grep { $edit_distances[$_] == $min_distance } 0..$#edit_distances;

#                 if (@min_distance_indices > 1) {
#                     $r_feature = $right{'line'} . q{:} . $rsequence;
#                 }
#                 else {
#                     $r_seqname   = $l_seqname;
#                     $r_feature   = $right{'line'} . q{:} . $rsequence;
#                     $r_start     = $reference{$l_seqname} + 1 - $left{'coord'} + $min_distance_indices[0];
#                     $r_end       = $reference{$l_seqname} + 1 - $left{'coord'} + $min_distance_indices[0] + $readsize;
#                     $r_score     = 1;
#                     $r_strand    = q{-};
#                     $r_attribute = 'target='
#                     . substr(
#                         $right_target_sequence,
#                         $min_distance_indices[0] - 2,
#                         $readsize + 4
#                     );
#                 }
#             }
        }
    }

    # multiple matches on left sequence
    if ( $lmatch > 1 && $rmatch == 1 ) {

        $l_source  = 'R/U';
        $r_source  = 'R/U';

        # gets the right chromosome
        $right{'chr0'} =~ m/(.*)/i;
        my $tmp = $1;
        $tmp =~ tr/A-Z/a-z/;
        $tmp =~ s/rc_//i;

        if ( $right{'chr0'} =~ m/^RC_/ ) {
            $r_seqname = $tmp;
            $r_feature = $right{'line'} . q{:} . $rsequence;
            $r_start   = $right{'coord'};
            $r_end     = $right{'coord'} + $rreadsize - 1;
            $r_score   = 1;
            $r_strand  = q{-};
            $r_frame   = $right{mm};
            $r_attribute = 'target='
                . substr(
                    $reference{"$tmp-rc"},
                    $right{'coord'} - 1,
                    $rreadsize + 4
                );
        }
        else {
            $r_seqname   = $tmp;
            $r_feature   = $right{'line'} . q{:} . $rsequence;
            $r_start     = $right{'coord'};
            $r_end       = $right{'coord'} + $rreadsize - 1;
            $r_score     = 1;
            $r_strand    = q{+};
            $r_frame     = $right{mm};
            $r_attribute = 'target='
                . substr(
                    $reference{"$tmp-seq"},
                    $right{'coord'} - 1,
                    $rreadsize + 4
                );
        }

        # loops through every possible left match
        my ( $bestcoord, $beststrand, $bestmm ) = (-1, -1, -1);
        my $score = 100000000; # arbitrary large number; lower is better

        U0_MM_LOOP_LEFT:
        for my $i ( 0 .. $lmatch - 1 ) {

            # checks best match
            my $tmpscore =
                check_mult_match ($left{"chr$i"}, $right{'chr0'}, $left{"coord$i"},
                                $right{'coord'}, $reference{$tmp}, $offset, $distance, $single_ends, $lreadsize, $rreadsize);

            # if $tmpscore is higher than current score of -1 discard potential match
            if ( $tmpscore == -1 or $tmpscore > $score ) {next U0_MM_LOOP_LEFT}

            # if $tmpscore flag is lower that current best score, keep
            if ( $tmpscore < $score ) {
                $score      = $tmpscore;
                $bestcoord  = $left{"coord$i"};
                $beststrand = $left{"chr$i"};
                $bestmm     = $left{"mm$i"};
            }
        }

        # gets the best left chromosome
        $beststrand =~ m/(.*)/i;
        $tmp = $1;
        $tmp =~ tr/A-Z/a-z/;
        $tmp =~ s/rc_//i;

        # if left sequence is distanced from right sequence within given range
        if ( $score >= $offset and $score <= $offset + $distance ) {

            if ( $beststrand =~ m/^RC_/i ) { # if left sequence maps to reverse strand

                $l_seqname = $tmp;
                $l_feature = $left{'line'} . q{:} . $lsequence;
                $l_start   = $bestcoord;
                $l_end     = $bestcoord + $lreadsize - 1;
                $l_score   = 1;
                $l_strand  = q{-};
                $l_frame   = $bestmm,
                $l_attribute = 'target='
                    . substr( $reference{"$tmp-rc"},
                              $bestcoord - 1,
                              $lreadsize + 4
                    );
            }
            else { # if left sequence maps to forward strand

                $l_seqname   = $tmp;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $bestcoord;
                $l_end       = $bestcoord + $lreadsize - 1;
                $l_score     = 1;
                $l_strand    = q{+};
                $l_frame     = $bestmm,
                $l_attribute = 'target='
                    . substr(
                        $reference{"$tmp-seq"},
                        $bestcoord - 1,
                        $lreadsize + 4
                    );
            }
        }
        else { # if right sequence is distanced from left sequence outside given range
            $l_feature = $left{'line'} . q{:} . $lsequence;
#             if ( $r_strand eq q{-} ) {

#                 my $left_target_sequence = substr $reference{"$r_seqname-rc"}, $right{'coord'} - 1 + $offset, $distance;
#                 $left_target_sequence =~ tr/Cc/Tt/;
#                 my @edit_distances = @{get_edit_distances ($lsequence, $left_target_sequence)};

#                 my $min_distance = min @edit_distances;
#                 my (@min_distance_indices) = grep { $edit_distances[$_] == $min_distance } 0..$#edit_distances;

#                 if (@min_distance_indices > 1) {
#                     $l_feature = $left{'line'} . q{:} . $lsequence;
#                 }
#                 else {
#                     $l_seqname   = $r_seqname;
#                     $l_feature   = $right{'line'} . q{:} . $lsequence;
#                     $l_start     = $right{'coord'} + $min_distance_indices[0];
#                     $l_end       = $right{'coord'} + $min_distance_indices[0] + $readsize;
#                     $l_score     = 1;
#                     $l_strand    = q{+};
#                     $l_attribute = 'target='
#                     . substr(
#                         $left_target_sequence,
#                         $min_distance_indices[0] - 2,
#                         $readsize + 4
#                     );
#                 }
#             }
#             else {

#                 my $left_target_sequence = substr $reference{"$r_seqname-seq"}, $right{'coord'} - 1 + $offset, $distance;
#                 $left_target_sequence =~ tr/Gg/Aa/;
#                 my @edit_distances = @{get_edit_distances ($lsequence, $left_target_sequence)};

#                 my $min_distance = min @edit_distances;
#                 my (@min_distance_indices) = grep { $edit_distances[$_] == $min_distance } 0..$#edit_distances;

#                 if (@min_distance_indices > 1) {
#                     $l_feature = $left{'line'} . q{:} . $lsequence;
#                 }
#                 else {
#                     $l_seqname   = $r_seqname;
#                     $l_feature   = $left{'line'} . q{:} . $lsequence;
#                     $l_start     = $reference{$r_seqname} + 1 - $right{'coord'} + $min_distance_indices[0];
#                     $l_end       = $reference{$r_seqname} + 1 - $right{'coord'} + $min_distance_indices[0] + $readsize;
#                     $l_score     = 1;
#                     $l_strand    = q{-};
#                     $l_attribute = 'target='
#                     . substr(
#                         $left_target_sequence,
#                         $min_distance_indices[0] - 2,
#                         $readsize + 4
#                     );
#                 }
#             }
        }
    }
    ##### One match on one end, multiple matches on other end #####




    ##### Multiple matches on both ends #####
    if ( $lmatch > 1 && $rmatch > 1 ) {

        $l_source  = 'R/R';
        $r_source  = 'R/R';

        # loops through each possible target on the left
        my ( $lbestcoord, $rbestcoord, $lbestmm, $lbeststrand, $rbeststrand, $rbestmm )
        = (-1, -1, -1, -1, -1, -1);
        my $bestscore = 100000000; # arbitrary large number; lower is better

        MM_LOOP_LEFT:
        for my $i ( 0 .. $lmatch - 1 ) {

            # gets current left chromosome
            $left{"chr$i"} =~ m/(.*)/i;
            my $tmp = $1;
            $tmp =~ tr/A-Z/a-z/;
            $tmp =~ s/rc_//i;

            # loops through each possible target on the right
            my ( $bestcoord, $beststrand, $bestmm ) = (-1, -1, -1);
            my $score = 100000000; # arbitrary large number; lower is better

            MM_LOOP_RIGHT:
            for my $j ( 0 .. $rmatch - 1 ) {

                my $tmpscore =
                    check_mult_match( $left{"chr$i"}, $right{"chr$j"},
                                    $left{"coord$i"}, $right{"coord$j"}, $reference{$tmp},
                                    $offset, $distance, $single_ends, $lreadsize, $rreadsize);

                # if $tmpscore if higher than current best score or invalid, discard current right sequence
                if ( $tmpscore == -1 or $tmpscore > $score ) {next MM_LOOP_RIGHT}

                # if $tmpscore is lower than current best score, keep current right sequence
                if ( $tmpscore < $score ) {
                    $score      = $tmpscore;
                    $bestcoord  = $right{"coord$j"};
                    $beststrand = $right{"chr$j"};
                    $bestmm     = $right{"mm$j"};
                }
            }

            # if no good matches on right side, skip this left sequence
            if ( $bestcoord == -1 && $beststrand == -1 ) {next MM_LOOP_LEFT}

            # if current score better than best score, keep current combination
            if ( $score < $bestscore ) {
                $bestscore   = $score;
                $lbestcoord  = $left{"coord$i"};
                $lbeststrand = $left{"chr$i"};
                $lbestmm     = $left{"mm$i"};
                $rbestcoord  = $bestcoord;
                $rbeststrand = $beststrand;
                $rbestmm     = $bestmm;
            }
        }

        # gets the current left chromosome
        $lbeststrand =~ m/(.*)/i;
        my $tmp = $1;
        $tmp =~ tr/A-Z/a-z/;
        $tmp =~ s/rc_//i;

        # if right sequence is distanced from left sequence within given range
        if ( $bestscore >= $offset and $bestscore <= $offset + $distance ) {

            # if the left sequence maps to the reverse strand
            if ( $lbeststrand =~ m/^RC_/i ) {

                $l_seqname   = $tmp;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $lbestcoord;
                $l_end       = $lbestcoord + $lreadsize - 1;
                $l_score     = 1;
                $l_strand    = q{-};
                $l_frame     = $lbestmm;
                $l_attribute = 'target='
                    . substr(
                        $reference{"$tmp-rc"},
                        $lbestcoord - 1,
                        $lreadsize + 4
                    );

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $rbestcoord;
                $r_end       = $rbestcoord + $rreadsize - 1;
                $r_score     = 1;
                $r_frame     = $rbestmm;

                if ($single_ends) {
                    $r_strand    = q{-};
                    $r_attribute = 'target='
                    . substr(
                        $reference{"$tmp-rc"},
                        $rbestcoord - 1,
                        $rreadsize + 4
                    );
                }
                else {
                    $r_strand    = q{+};
                    $r_attribute = 'target='
                    . substr(
                        $reference{"$tmp-seq"},
                        $rbestcoord - 1,
                        $rreadsize + 4
                    );
                }
            }
            else { # if the left sequence maps to the reverse strand

                $l_seqname   = $tmp;
                $l_feature   = $left{'line'} . q{:} . $lsequence;
                $l_start     = $lbestcoord;
                $l_end       = $lbestcoord + $lreadsize - 1;
                $l_score     = 1;
                $l_strand    = q{+};
                $l_frame     = $lbestmm;
                $l_attribute = 'target='
                    . substr(
                        $reference{"$tmp-seq"},
                        $lbestcoord - 1,
                        $lreadsize + 4
                    );

                $r_seqname   = $tmp;
                $r_feature   = $right{'line'} . q{:} . $rsequence;
                $r_start     = $rbestcoord;
                $r_end       = $rbestcoord + $rreadsize - 1;
                $r_score     = 1;
                $r_frame     = $rbestmm;

                if ($single_ends) {
                    $r_strand    = q{+};
                    $r_attribute = 'target='
                    . substr(
                        $reference{"$tmp-seq"},
                        $rbestcoord - 1,
                        $rreadsize + 4
                    );
                }
                else {
                    $r_strand    = q{-};
                    $r_attribute = 'target='
                    . substr(
                        $reference{"$tmp-rc"},
                        $rbestcoord - 1,
                        $rreadsize + 4
                    );
                }
            }
        }
        else { # if right sequence is distanced from left sequence outside given range
            $l_source = "R/R";
            $r_source = "R/R";

            $l_feature = $left{'line'} . q{:} . $lsequence;
            $r_feature = $right{'line'} . q{:} . $rsequence;
        }
    }
    ##### Multiple matches on both ends #####


    ##### END POSSIBLE CASES HERE #####


    $r_strand = $l_strand if $single_ends;

    # print both ends
    print join( "\t",
                $l_seqname, $l_source, $l_feature, $l_start, $l_end,
                $l_score,   $l_strand, $l_frame,  $l_attribute ),
                "\n"
                unless $skip_nm and $l_score = 1;
    print join( "\t",
                $r_seqname, $r_source, $r_feature, $r_start, $r_end,
                $r_score,   $r_strand, $r_frame,  $r_attribute ),
                "\n"
                unless $skip_nm and $r_score = 1;
}

# end for loop through all sequences
print STDERR "done!\n";
close STDOUT;
close $LEFT;
close $RIGHT;

# main program finished




##### Start function definitions #####

# Given 2 chromosomes and 2 coordinates, plus a library offset and distance (see header for details):
# Checks that each element of a pair belongs to the reverse complement of the other,
# and that both their coordinates summed don't exceed the length of the chromosome plus or minus the library size.
# Returns: $match = 1 or 0
sub check_match {
    my ( $lchr, $rchr, $lcoord, $rcoord, $chrlength, $offset, $distance, $single_ends, $lreadsize, $rreadsize ) = @_;

    if ($single_ends) {
        # check chromosome
        if ( $lchr =~ m/^([^R][^C][^_])/i and $rchr =~ m/^RC_$1/i) {
            return 0;
        }
        elsif ( $rchr =~ m/^RC_(.*)/i and $lchr =~ m/^$1/i ) {
            return 0;
        }
    }
    else {
        if ( $lchr =~ m/^([^R][^C][^_])/i and $rchr !~ m/^RC_$1/i) {
            return 0;
        }
        elsif ( $rchr =~ m/^RC_(.*)/i and $lchr !~ m/^$1/i ) {
            return 0;
        }
    }

    # check coordinates
    my $coord_distance = abs (($chrlength - $lcoord + 1) - $rcoord);
    
    $coord_distance = abs ($lcoord - $rcoord) - $lreadsize if $single_ends;

    if ( $coord_distance > ($offset + $distance) ) {return 0;}

    return 1;
}

# Given 2 chromosomes and 2 coordinates, plus a library offset and distance (see header for details):
# Checks that each element of a pair belongs to the reverse complement of the other,
# and that both their coordinates summed don't exceed the length of the chromosome plus or minus the library size.
# Returns $score = library size/distance between both sequences or -1 for out of range distances
sub check_mult_match {
    my ( $lchr, $rchr, $lcoord, $rcoord, $chrlength, $offset, $distance, $single_ends, $lreadsize, $rreadsize ) = @_;
    my $score;

    if ($single_ends) {
        # check chromosome
        if ( $lchr =~ m/^([^R][^C][^_])/i and $rchr =~ m/^RC_$1/i) {
            return -1;
        }
        elsif ( $rchr =~ m/^RC_(.*)/i and $lchr =~ m/^$1/i ) {
            return -1;
        }
    }
    else {
        if ( $lchr =~ m/^([^R][^C][^_])/i and $rchr !~ m/^RC_$1/i) {
            return -1;
        }
        elsif ( $rchr =~ m/^RC_(.*)/i and $lchr !~ m/^$1/i ) {
            return -1;
        }
    }

    my $coord_distance = abs (($chrlength - $lcoord + 1) - $rcoord);

    $coord_distance = abs ($lcoord - $rcoord) - $lreadsize if $single_ends;

    if ( $coord_distance > ($offset + $distance) ) {return -1;}

    return $coord_distance;
}

# Parses single eland3 record (output from seqmap)
# Returns reference to hash with following keys:
# 'line', 'sequence', 'matches', 'chr[0-#matches-1]', 'coord[0-#matches-1]'
sub parseEland3Line {
    my $tmp_eland_line = shift;
    $tmp_eland_line =~ s/[\n\r]//g;
    my @eland_line = (split /\t/, $tmp_eland_line);

    my %hash = ();
    $hash{'line'}     = $eland_line[0];
    $hash{'sequence'} = $eland_line[1];

    if ($eland_line[2] =~ m/NM/) {
        $hash{'matches'} = 0;
    }

    elsif ($eland_line[2] =~ m/^[0-9]+$/) {
        $hash{'matches'} = $eland_line[2];
    }

    elsif ($eland_line[2] =~ m/:/) {
        my @all_matches = split /,/, $eland_line[3];
        $hash{'matches'} = scalar @all_matches;
    }

    if ($hash{'matches'} > 1) {
        my @all_reads = split /,/, $eland_line[3];

        for my $i (0..@all_reads - 1) {
            ($hash{"chr$i"}, $hash{"coord$i"}) = split /:/, $all_reads[$i];
            $hash{"coord$i"} =~ s/[A-Z]([0-9])$//i;
            $hash{"mm$i"} = $1;
        }
    }
    elsif ($hash{'matches'} == 1) {
        ($hash{'chr0'}, $hash{'coord'}) = split /:/, $eland_line[3];
        $hash{'coord'} =~ s/[A-Z]([0-9])$//i;
        $hash{'mm'} = $1;
    }
    return \%hash;
}

# Converts input sequence to reverse complement
# Returns scalar $string with processed sequence
sub reverseComp {
    my $shortseq = shift;
    $shortseq =~ tr/ACGTacgt/TGCAtgca/;
    $shortseq =~ s/\n//;
    return reverse $shortseq;
}


sub generate_windows {
    my ($target_sequence, $readsize) = @_;

    my @targets = ();

    for (my $i = 0; $i < (length $target_sequence) - $readsize; $i++) {
        push @targets, substr $target_sequence, $i, $readsize;
    }

    return \@targets;
}


sub get_edit_distances {
    my ($base_sequence, $target_sequence) = @_;

    my $targets_ref = generate_windows ($target_sequence, length $base_sequence);

    my @edit_distances;

    for my $target (@{$targets_ref}) {
        push @edit_distances, distance ($base_sequence, $target);
    }

    return \@edit_distances;
}


sub usage {
    if ( $usage || @ARGV < 5 ) {
        print STDERR "correlatePairedEnds.pl <PARAMETERS> [OPTIONS]
\t<--left>\t5' raw sequence file
\t<--right>\t3' raw sequence file
\t<--reference>\tFasta genome file
\t<--offset>\tMinimum library size
\t<--distance>\tMaximum variation from offset
\t[--output]\tFilename to write results to (default is STDOUT)
\t[--readsize]\tRaw sequence length
\t[--single-ends]
\t[--quiet]\tDon't output perl warnings
\t[--verbose]\tOutput detailed perl diagnostic messages
\t[--usage]\tPrints this
Takes a paired ends sequencing output (needs to be preprocessed to fit above) 
and compares each pair to maximize unique alignment matches to a genome.\n";
        exit 1;
    }
    return 0;
}
