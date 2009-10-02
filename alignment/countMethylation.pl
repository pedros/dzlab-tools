#!/usr/bin/perl
#
# Last edited 2009-08-25
# Copyright 2008-2009 Pedro Silva <psilva@nature.berkeley.edu/>
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

Takes a paired ends alignment output in gff form (from correlatePairedEnds.pl) and 
outputs a gff file in which each record is an individual cytosin bp in the scaffold sequence.

=cut

use Data::Dumper;
use Getopt::Long;
use strict;
use warnings;
use diagnostics;
disable diagnostics;
use List::Util qw(max sum);
use Carp;
#use Smart::Comments '###';

# Globals, passed as command line options
my $gfffile = '';
my $reference = '';
my $stats_only = 0;
my $sort = 0;
my $di_nucleotide_count = 0;
my $output = "-";
my $verbose = 0;
my $quiet = 0;
my $usage = 0;

# Initial check of command line parameters
usage ();

# Grabs and parses command line options
my $result = GetOptions (
    'gff|f=s'      => \$gfffile,
    'ref|r=s'      => \$reference,
    'sort|s'       => \$sort,
    'di-nuc|d:o'   => \$di_nucleotide_count,
    "output|o:s"   => \$output,
    'stats-only|t' => \$stats_only,
    "verbose|v"    => sub {enable diagnostics;},
    "quiet|q"      => sub {disable diagnostics;no warnings;},
    "usage|help|h" => \&usage
    );

$output = '/dev/null' if $stats_only;

# redirects STDOUT to file if specified by user
if(!($output eq '-')) {
    open(STDOUT, '>', "$output") or die("Can't redirect STDOUT to file: $output");
}

if ($sort) {
    # opens gff file
    open my $GFFIN, '<', $gfffile or die("Can't read file: $gfffile");
    my @gff_data = <$GFFIN>;
    close $GFFIN;

    open my $GFFOUT, '>', $gfffile or die("Can't write to file: $gfffile");
    @gff_data =  gff_sort (\@gff_data);
    print $GFFOUT @gff_data;
    close $GFFOUT;
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
	   "'C' count; 'T' count",
    ), "\n";

# declare and initialize a hash of hashes, where each inner hash's key is a 'C' coordinate
# and its key/value pairs are total c/t counts, and contexts
my %HoH=();

# for basic stats
my %total_count = ();
($total_count{bp}, $total_count{overlaps}) = (0, 0);

# holds name of chromosomes as keys and length of chromosomes in bp as values
my %reference = ();

# begins limited scope for extracting names and lengths of chromosomes in reference file
{
    # reads in the reference genome file into @fastaseq
    open my $REF, '<', "$reference" or croak "Can't open file: $reference";
    my @fastaseq = <$REF>;
    close $REF;

    # find and store indices for each chromosome change and corresponding descriptions
    my ( @idx, @dsc ) = ();
    for my $i ( 0 .. @fastaseq - 1 ) { ### Indexing $reference...  % done
        if ( $fastaseq[$i] =~ m/^>/ ) {
            $fastaseq[$i] =~ s/>//g;
            $fastaseq[$i] = ( split /\s/, "$fastaseq[$i]" )[0];
            $fastaseq[$i] =~ tr/A-Z/a-z/;
            push @idx, $i;
            push @dsc, $fastaseq[$i];
        }
    }

    # gets and saves each chromosome's sequence and reverse complemented sequence
    for my $j ( 0 .. @idx - 1 ) { ### Loading $reference into memory...  % done
        my $line;
        if ( $j == scalar @idx - 1 ) {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. @fastaseq - 1]);
        }
        else {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. $idx[$j + 1] - 1]);
        }
        $line =~ s/[\n\r]//g;
        $reference{ $dsc[$j] } = length $line;
    }
}


# opens gff file
open my $GFF, "<", $gfffile or die("Can't read file: $gfffile");

my $last_seen_coordinate = 1;

my %reads_buffer = ();

# reads in input file until EOF
GFF_READ_LOOP:
while (<$GFF>) { ### Indexing...
    chomp $_;

    # skips comments ('#') and blank lines
    next GFF_READ_LOOP if ($_ =~ m/^#.*$|^\s*$/);

    # parses each record line into a simple hash
    my %record = %{ readGFF ($_) };

    # skips records with non-unique matches
    next GFF_READ_LOOP if ($record{'score'} != 1);

     if (scalar keys %HoH && $record{'start'} > $last_seen_coordinate) {
         sort_and_count (\%HoH, $di_nucleotide_count);
     }

    $last_seen_coordinate = $record{'end'};

    # grabs methylated sequence from gff file
    # assumes sequence is exactly $readsize bps long
    $record{'feature'} =~ m/([ACGTN]+$)/ or die Dumper \%record;
    my @methylated = split(//, $1);

    # keep track of total bp seen
    $total_count{bp} += scalar @methylated;

    # grabs scaffold sequence from gff file
    # assumes sequence is in the 'attribute' field, separated from extraneous information by a '='
    my @unmethylated = split(//, (split "=", $record{'attribute'})[1]);

    # check for overlapping pairs of reads
    my ($pair_id, $read_id) = $record{'feature'} =~ m/(^.*?)\/([12]):[ACGTN]+/;
    my $overlap = 0; # array reference with overlap start coord, end coord

    # if matching pair has not been seen, create a new entry for this pair in buffer
    if (!exists $reads_buffer{$pair_id}) {
        $reads_buffer{$pair_id} = [$record{'start'}, $record{'end'}, $read_id];
    }
    else { # if matching pair on buffer
        # existing end runs into pair - overlap is current pair's start to existing end's end
        if ($record{'end'} > $reads_buffer{$pair_id}->[0] and $record{'start'} < $reads_buffer{$pair_id}->[1]) {
            $overlap = [$record{'start'}, $reads_buffer{$pair_id}->[1]];
        }
        # current end runs into pair - overlap is existing end's start to current end's end
        elsif ($reads_buffer{$pair_id}->[1] > $record{'start'} and $reads_buffer{$pair_id}->[0] < $record{'end'}) {
            $overlap = [$reads_buffer{$pair_id}->[0], $record{'end'}];
        }
        delete $reads_buffer{$pair_id};
    }

    # loops through each character in current sequences
    READ:
    for ( my $i = $record{'start'}; $i < $record{'end'}; $i++ ) {

	# sets $j to 0 for forward strand coordinates
	my $j = $i - $record{'start'};

	# Since the input coordinates are absolute (ie. from the forward strand)
	# if a particular read maps to the reverse strand, the coordinates for any
	# 'c's hit need to be reversed
	my $coord = $i;
	$coord = $reference{$record{'seqname'}} - $i + 1 if $record{'strand'} eq '-';

        # checks that we're looking at a left/1 sequence AND current character is a 'C'
	# this regex is adapted to the Solexa sequences we have, so it might change in the future
	# we're looking 2 characters ahead in the scaffold because it has 2 extra chars in the beginning
 	if ( $unmethylated[$j + 2] =~ m/[Cc]/
             and $record{'feature'} =~ m/\/1:[ACGTN]+/) {

	    # checks what happened in the methylated sequence
	    # and updates the appropriate 'c' or 't' count
	    if ( $methylated[$j] =~ m/[Cc]/ ) {
		$HoH{$coord}[0]++; # c_count
	    }
	    elsif ($methylated[$j] =~ m/[Tt]/ ) {
		$HoH{$coord}[1]++; # t_count
	    }

	    # checks the context by looking ahead +3 or +3..+4 bps
	    # because the scaffold read is displaced by 2 bps
	    if ( $unmethylated[$j+3] =~ m/[Gg]/ ) {
                $HoH{$coord}[2]++; # cg_count
	    }
	    elsif ( join("", @unmethylated[ ($j + 3)..($j + 4) ] ) =~ m/[^Gg][Gg]/) {
		$HoH{$coord}[3]++; # chg_count
	    }
	    elsif ( join("", @unmethylated[ ($j + 3)..($j + 4) ] ) =~ m/[^Gg][^Gg]/) {
		$HoH{$coord}[4]++; # chh_count
	    }

            if ($di_nucleotide_count) {
                # count di-nucleotide contexts
                if ( $unmethylated[$j+3] =~ m/[Aa]/ ) {
                    $HoH{$coord}[11]++; # ca_count
                }
                elsif ( $unmethylated[$j+3] =~ m/[Cc]/ ) {
                    $HoH{$coord}[12]++; # cc_count
                }
                elsif ( $unmethylated[$j+3] =~ m/[Tt]/ ) {
                    $HoH{$coord}[13]++; # ct_count
                }
            }

	    # grab some necessary information into the data structure
	    # we will print this information later
	    $HoH{$coord}[5] = $coord; # coord
	    $HoH{$coord}[6] = $record{'seqname'}; # chromosome
	    $HoH{$coord}[7] = $record{'strand'}; # strand
            $HoH{$coord}[9]++ if $record{'strand'} eq q{+};
            $HoH{$coord}[10]++ if $record{'strand'} eq q{-};
	}

	# checks that we're looking at a right/2 sequence AND current character is a 'C'
	# this regex is adapted to the Solexa sequences we have, so it might change in the future
	# we're looking 2 characters before in the scaffold because it has 2 extra chars in the beginning
	# AND we're looking in the reverse strand, so the context goes right-to-left
	elsif ( $unmethylated[$j + 2] =~ m/[Gg]/
                and $record{'feature'} =~ m/\/2:[ACGTN]+/) {

            # checks what happened in the methylated sequence
	    # and updates the appropriate 'c' or 't' count
	    if ( $methylated[$j] =~ m/[Gg]/ ) {
		$HoH{$coord}[0]++; # c_count
	    }
	    elsif ( $methylated[$j] =~ m/[Aa]/ ) {
		$HoH{$coord}[1]++; # t_count
	    }

	    # checks the context by looking behind +1 or +0..+1 bps
	    # because the scaffold read is displaced by 2 bps
	    if ( $unmethylated[$j + 1] =~ m/[Cc]/ ) {
                $HoH{$coord}[2]++; # cg_count
	    }
	    elsif ( join("", @unmethylated[ $j..($j + 1) ] ) =~ m/[Cc][^Cc]/ ) {
		$HoH{$coord}[3]++; # chg_count
	    }
	    elsif ( join("", @unmethylated[ $j..($j + 1) ] ) =~ m/[^Cc][^Cc]/ ) {
		$HoH{$coord}[4]++; # chh_count
	    }

            if ($di_nucleotide_count) {
                # count di-nucleotide contexts
                if ( $unmethylated[$j + 1] =~ m/[Tt]/ ) {
                    $HoH{$coord}[11]++; # ca_count
                }
                elsif ( $unmethylated[$j + 1] =~ m/[Gg]/ ) {
                    $HoH{$coord}[12]++; # cc_count
                }
                elsif ( $unmethylated[$j + 1] =~ m/[Aa]/ ) {
                    $HoH{$coord}[13]++; # ct_count
                }
            }

            # if current coordinate is in pair overlaping region only count /1 reads
            if ($overlap and
                $overlap->[0] < $coord and
                $overlap->[1] > $coord) {
                next GFF_READ_LOOP;
                $total_count{overlaps}++;
            }

	    # grab some necessary information into the data structure
	    # we will print this information later
	    $HoH{$coord}[5] = $coord; # coord
	    $HoH{$coord}[6] = $record{'seqname'}; # chromosome
	    $HoH{$coord}[7] = q{+} if $record{'strand'} eq q{-}; # strand
	    $HoH{$coord}[7] = q{-} if $record{'strand'} eq q{+}; # strand
            $HoH{$coord}[9]++ if $HoH{$coord}[7] eq q{+};
            $HoH{$coord}[10]++ if $HoH{$coord}[7] eq q{-};
	}
    }
}

sort_and_count (\%HoH, $di_nucleotide_count) if scalar keys %HoH > 0;
close($GFF);

# open for counting frequencies
my @freq = count_freq($output, $di_nucleotide_count);

open my $FREQ_OUT, '>', "$output.freq" or die "Can't write to file: $output.freq";

unless ($di_nucleotide_count) {
    print $FREQ_OUT join ("\t",
                          'bp',
                          'overlaps',
                          'C',
                          'CG',
                          'CHG',
                          'CHH',
                          'T',
                          'TG',
                          'THG',
                          'THH',
                          'C_ratio',
                          'CG_ratio',
                          'CHG_ratio',
                          'CHH_ratio',
                          'filtered_C',
                          'filtered_CG',
                          'filtered_CHG',
                          'filtered_CHH',
                          'filtered_T',
                          'filtered_TG',
                          'filtered_THG',
                          'filtered_THH',
                          'filtered_C ratio',
                          'filtered_CG ratio',
                          'filtered_CHG ratio',
                          'filtered_CHH ratio',
                      ), "\n";
}
else {
    print $FREQ_OUT join ("\t",
                          'bp',
                          'overlaps',
                          'C',
                          'NO-CG',
                          'CG',
                          'CA',
                          'CC',
                          'CT',
                          'T',
                          'NO-TG',
                          'TG',
                          'TA',
                          'TC',
                          'TT',
                          'C_ratio',
                          'CG_ratio',
                          'CA_ratio',
                          'CC_ratio',
                          'CT_ratio',
                          'filtered_C',
                          'filtered-NO-CG',
                          'filtered_CG',
                          'filtered_CA',
                          'filtered_CC',
                          'filtered_CT',
                          'filtered_T',
                          'filtered-NO-TG',
                          'filtered_TG',
                          'filtered_TA',
                          'filtered_TC',
                          'filtered_TT',
                          'filtered_C ratio',
                          'filtered_CG ratio',
                          'filtered_CA ratio',
                          'filtered_CC ratio',
                          'filtered_CT ratio',
                      ), "\n";
}


print $FREQ_OUT join ("\t",
                      $total_count{bp},
                      $total_count{overlaps},
                      @freq,
                  ), "\n";
close $FREQ_OUT;

# done



sub gff_sort {
    my $data_ref = shift;
    return sort {
        (split /\t/, $a)[3] <=> (split /\t/, $b)[3]
    } @{$data_ref};
}


sub sort_and_count {
    my ($single_count_ref, $di_nucleotide_count)
    = @_;
    my %single_count = %{$single_count_ref};

    # loops through every initialized key in main hash
    # keys are sorted, so the output is going to be sorted by starting coordinate
    for my $i (sort {$a <=> $b} keys %single_count) { ### Sorting and counting...  % done

        # we need to check that a given key/value pair was initialized
        # if it wasn't, initialize it to zero (to avoid division-by-zero, etc)
        for (0..4) {
            $single_count{$i}[$_] = 0 if !exists $single_count{$i}[$_];
        }

        # calculates c/(c+t), which varies between 0 and 1
        # 0 is no methylation, 1 is total methylation
        if ($single_count{$i}[0] + $single_count{$i}[1] == 0) {
            $single_count{$i}[8] = 0;
        }
        else {
            $single_count{$i}[8] = $single_count{$i}[0]/($single_count{$i}[0]+$single_count{$i}[1]);
        }

        # finds the appropriate context to put in 'feature' field
        my $context = q{.};

        if ($di_nucleotide_count) {
            # we need to check that a given key/value pair was initialized
            # if it wasn't, initialize it to zero (to avoid division-by-zero, etc)
            for (11 .. 13) {
                $single_count{$i}[$_] = 0 if !exists $single_count{$i}[$_];
            }

            if ( $single_count{$i}[2] > 0) {
                if ($single_count{$i}[2] > $single_count{$i}[11] && $single_count{$i}[2] > $single_count{$i}[12] && $single_count{$i}[2] > $single_count{$i}[13]) {
                    $context = 'CG';
                }
            }
            elsif ($single_count{$i}[11] > 0) {
                if ($single_count{$i}[11] > $single_count{$i}[12] && $single_count{$i}[11] > $single_count{$i}[13]) {
                    $context = 'CA';
                }
            }
            elsif ($single_count{$i}[12] > 0) {
                if ($single_count{$i}[12] > $single_count{$i}[13]) {
                    $context = 'CC';
                }
            }
            elsif ($single_count{$i}[13]>0) {
                $context = 'CT';
            }

            # checks that the context for any given 'c' is determinable
            # if not just output a '.' in its place
            if (($single_count{$i}[2] > 0 && $single_count{$i}[11] > 0) ||
                ($single_count{$i}[2] > 0 && $single_count{$i}[12] > 0) || 
                ($single_count{$i}[2] > 0 && $single_count{$i}[13] > 0) ||
                ($single_count{$i}[11] > 0 && $single_count{$i}[12] > 0) ||
                ($single_count{$i}[11] > 0 && $single_count{$i}[13] > 0) ||
                ($single_count{$i}[12] > 0 && $single_count{$i}[13] > 0)) {
                $context = q{.};
            }
        }
        else {
            if ( $single_count{$i}[2] > 0) {
                if ($single_count{$i}[2]>$single_count{$i}[3] && $single_count{$i}[2]>$single_count{$i}[4]) {
                    $context = 'CG';
                }
            }
            elsif ($single_count{$i}[4] > 0) {
                if ($single_count{$i}[4] > $single_count{$i}[3]) {
                    $context = 'CHH';
                }
            }
            elsif ($single_count{$i}[3]>0) {
                $context = 'CHG';
            }

            # checks that the context for any given 'c' is determinable
            # if not just output a '.' in its place
            if (($single_count{$i}[2] > 0 && $single_count{$i}[3] > 0) ||
                ($single_count{$i}[2] > 0 && $single_count{$i}[4] > 0) ||
                ($single_count{$i}[3] > 0 && $single_count{$i}[4] > 0)) {
                $context = q{.};
            }
        }


        if (defined $single_count{$i}[9] and defined $single_count{$i}[10]) {
            $single_count{$i}[7] = q{.};
        }

        my $attribute = join q{;}, "c=$single_count{$i}[0]", "t=$single_count{$i}[1]";
        $attribute .= '*' if ($single_count{$i}[0] + $single_count{$i}[1] > 20
        and $single_count{$i}[6] !~ m/chr[cm]/);

        # prints a single gff record
        print join("\t",
                   $single_count{$i}[6],
                   'dz_cm',
                   $context,
                   $single_count{$i}[5],
                   $single_count{$i}[5],
                   sprintf("%e", $single_count{$i}[8]), # score
                   $single_count{$i}[7],
                   q{.},
                   $attribute,
               ), "\n";

        delete $single_count{$i};
    }
    %{$single_count_ref} = ();
}


# readGFF reads in a single GFF line, and outputs a hash reference
sub readGFF {
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute)
    = split(/\t/, $_[0]);

    my %rec = (
	'seqname'  =>$seqname,
	'source'   =>$source,
	'feature'  =>$feature,
	'start'    =>$start,
	'end'      =>$end,
	'score'    =>$score,
	'strand'   =>$strand,
	'frame'    =>$strand,
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
\t[--sort]\tPre-sort input file in-place before processing
\t[--output]\tFilename to write results to (default is STDOUT)
\t[--di-nuc]\tCount C di-nucleotide contexts instead of tri-nuc
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


sub count_freq {
    my ($gfffile, $di_nucleotide_count) = @_;
    open my $GFF, "<", $gfffile or die("Can't read file: $gfffile");

    my (%filtered, %unfiltered) = ();
    while (<$GFF>) {
        chomp;

        next if ($_ =~ m/^#.*$|^\s*$/);
        my ($context, $attribute) = (split /\t/, $_)[2, -1];
        my ($c, $t) = (split /;/, $attribute);

        my $filter = 1 if $t =~ m/\*$/;

        ($c) = $c =~ m/c=(\d+)/;
        ($t) = $t =~ m/t=(\d+)/;

        $unfiltered{$context}{c} = 0 unless defined $unfiltered{$context}{c};
        $unfiltered{$context}{t} = 0 unless defined $unfiltered{$context}{t};
        $filtered{$context}{c}   = 0 unless defined $filtered{$context}{c};
        $filtered{$context}{t}   = 0 unless defined $filtered{$context}{t};

        $unfiltered{$context}{c} += $c;
        $unfiltered{$context}{t} += $t;

        $filtered{$context}{c} += $c unless $filter;
        $filtered{$context}{t} += $t unless $filter;
    }
    close $GFF;

    return unless keys %filtered and keys %unfiltered;

    my ($total_unfiltered_c, $total_unfiltered_t, $total_filtered_c, $total_filtered_t);
    if ($di_nucleotide_count) {
        $total_unfiltered_c = $unfiltered{CG}{c} + $unfiltered{CA}{c} + $unfiltered{CC}{c} + $unfiltered{CT}{c};
        $total_unfiltered_t = $unfiltered{CG}{t} + $unfiltered{CA}{t} + $unfiltered{CT}{t} + $unfiltered{CT}{t};

        $total_filtered_c = $filtered{CG}{c} + $filtered{CA}{c} + $filtered{CC}{c} + $filtered{CT}{c};
        $total_filtered_t = $filtered{CG}{t} + $filtered{CA}{t} + $filtered{CC}{t} + $filtered{CT}{t};
    }
    else {
        $total_unfiltered_c = $unfiltered{CG}{c} + $unfiltered{CHG}{c} + $unfiltered{CHH}{c};
        $total_unfiltered_t = $unfiltered{CG}{t} + $unfiltered{CHG}{t} + $unfiltered{CHH}{t};
        
        $total_filtered_c = $filtered{CG}{c} + $filtered{CHG}{c} + $filtered{CHH}{c};
        $total_filtered_t = $filtered{CG}{t} + $filtered{CHG}{t} + $filtered{CHH}{t};
    }

    my (
        $unfiltered_C_ratio, $filtered_C_ratio,
        $unfiltered_CG_ratio, $filtered_CG_ratio,
        $unfiltered_CA_ratio, $filtered_CA_ratio,
        $unfiltered_CC_ratio, $filtered_CC_ratio,
        $unfiltered_CT_ratio, $filtered_CT_ratio,
        $unfiltered_CHG_ratio, $filtered_CHG_ratio,
        $unfiltered_CHH_ratio, $filtered_CHH_ratio,
    ) = 0 x 14;

    $unfiltered_C_ratio = $total_unfiltered_c / ($total_unfiltered_c + $total_unfiltered_t) if $total_unfiltered_c + $total_unfiltered_t != 0;
    $filtered_C_ratio = $total_filtered_c / ($total_filtered_c + $total_filtered_t) if $total_filtered_c + $total_filtered_t != 0;

    $unfiltered_CG_ratio = $unfiltered{CG}{c} / ($unfiltered{CG}{c} + $unfiltered{CG}{t}) if $unfiltered{CG}{c} + $unfiltered{CG}{t} != 0;
    $filtered_CG_ratio = $filtered{CG}{c} / ($filtered{CG}{c} + $filtered{CG}{t}) if $filtered{CG}{c} + $filtered{CG}{t} != 0;

    if ($di_nucleotide_count) {
        $unfiltered_CA_ratio = $unfiltered{CA}{c} / ($unfiltered{CA}{c} + $unfiltered{CA}{t}) if $unfiltered{CA}{c} + $unfiltered{CA}{t} != 0;
        $filtered_CA_ratio   = $filtered{CA}{c}   / ($filtered{CA}{c}   + $filtered{CA}{t}) if $filtered{CA}{c} + $filtered{CA}{t} != 0;

        $unfiltered_CC_ratio = $unfiltered{CC}{c} / ($unfiltered{CC}{c} + $unfiltered{CC}{t}) if $unfiltered{CC}{c} + $unfiltered{CC}{t} != 0;
        $filtered_CC_ratio   = $filtered{CC}{c}   / ($filtered{CC}{c}   + $filtered{CC}{t}) if $filtered{CC}{c} + $filtered{CC}{t} != 0;

        $unfiltered_CT_ratio = $unfiltered{CT}{c} / ($unfiltered{CT}{c} + $unfiltered{CT}{t}) if $unfiltered{CT}{c} + $unfiltered{CT}{t} != 0;
        $filtered_CT_ratio   = $filtered{CT}{c}   / ($filtered{CT}{c}   + $filtered{CT}{t}) if $filtered{CT}{c} + $filtered{CT}{t} != 0;
        
    }
    else {
        $unfiltered_CHG_ratio = $unfiltered{CHG}{c} / ($unfiltered{CHG}{c} + $unfiltered{CHG}{t}) if $unfiltered{CHG}{c} + $unfiltered{CHG}{t} != 0;
        $filtered_CHG_ratio   = $filtered{CHG}{c}   / ($filtered{CHG}{c}   + $filtered{CHG}{t}) if $filtered{CHG}{c} + $filtered{CHG}{t} != 0;

        $unfiltered_CHH_ratio = $unfiltered{CHH}{c} / ($unfiltered{CHH}{c} + $unfiltered{CHH}{t}) if $unfiltered{CHH}{c} + $unfiltered{CHH}{t} != 0;
        $filtered_CHH_ratio   = $filtered{CHH}{c}   / ($filtered{CHH}{c}   + $filtered{CHH}{t}) if $filtered{CHH}{c} + $filtered{CHH}{t} != 0;
    }


    if ($di_nucleotide_count) {
        return (
            $total_unfiltered_c,
            $total_unfiltered_c - $unfiltered{CG}{c},,
            $unfiltered{CG}{c},
            $unfiltered{CA}{c},
            $unfiltered{CC}{c},
            $unfiltered{CT}{c},
            $total_unfiltered_t,
            $total_unfiltered_t - $unfiltered{CG}{t},
            $unfiltered{CG}{t},
            $unfiltered{CA}{t},
            $unfiltered{CC}{t},
            $unfiltered{CT}{t},
            $unfiltered_C_ratio,
            $unfiltered_CG_ratio,
            $unfiltered_CA_ratio,
            $unfiltered_CC_ratio,
            $unfiltered_CT_ratio,
            $total_filtered_c,
            $total_filtered_c - $filtered{CG}{c},
            $filtered{CG}{c},
            $filtered{CA}{c},
            $filtered{CC}{c},
            $filtered{CT}{c},
            $total_filtered_t,
            $total_filtered_t - $filtered{CG}{t},
            $filtered{CG}{t},
            $filtered{CA}{t},
            $filtered{CC}{t},
            $filtered{CT}{t},
            $filtered_C_ratio,
            $filtered_CG_ratio,
            $filtered_CA_ratio,
            $filtered_CC_ratio,
            $filtered_CT_ratio,
        );
    }
    else {
        return (
            $total_unfiltered_c,
            $unfiltered{CG}{c},
            $unfiltered{CHG}{c},
            $unfiltered{CHH}{c},
            $total_unfiltered_t,
            $unfiltered{CG}{t},
            $unfiltered{CHG}{t},
            $unfiltered{CHH}{t},
            $unfiltered_C_ratio,
            $unfiltered_CG_ratio,
            $unfiltered_CHG_ratio,
            $unfiltered_CHH_ratio,
            $total_filtered_c,
            $filtered{CG}{c},
            $filtered{CHG}{c},
            $filtered{CHH}{c},
            $total_filtered_t,
            $filtered{CG}{t},
            $filtered{CHG}{t},
            $filtered{CHH}{t},
            $filtered_C_ratio,
            $filtered_CG_ratio,
            $filtered_CHG_ratio,
            $filtered_CHH_ratio,
        );
    }
}
