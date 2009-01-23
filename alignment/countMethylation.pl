#!/usr/bin/perl
#
# Last edited 2008-12-02
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

Takes a paired ends alignment output in gff form (from correlatePairedEnds.pl) and 
outputs a gff file in which each record is an individual cytosin bp in the scaffold sequence.

=cut

use Data::Dumper;
use Getopt::Long;
use strict;
use warnings;
use diagnostics;
disable diagnostics;
use List::Util qw(max);

# Globals, passed as command line options
my $gfffile = "";
my $sort = 0;
my $output = "-";
my $verbose = 0;
my $quiet = 0;
my $usage = 0;

# Initial check of command line parameters
usage ();

# Grabs and parses command line options
my $result = GetOptions (
    "gff|f=s" => \$gfffile,
    'sort|s'  => \$sort,
    "output|o:s" => \$output,
    "verbose|v" => sub {enable diagnostics;use Smart::Comments},
    "quiet|q" => sub {disable diagnostics;no warnings;},
    "usage|help|h" => \&usage
    );

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
# this becomes BIG (=~ 8GB for 400MB input files) and needs to be updated somehow.
# most likely the data structure will need to be swaped for a database (sqlite should be fine)
my %HoH=();

my %total_count = ();

($total_count{c_count},
 $total_count{t_count},
 $total_count{cg_count},
 $total_count{chg_count},
 $total_count{chh_count},
) = (0, 0, 0, 0, 0);
   
# opens gff file
open my $GFF, "<", $gfffile or die("Can't read file: $gfffile");

# reads in input file until EOF
GFF_READ_LOOP:
while(<$GFF>) { ### Indexing...
    chomp $_;

    # skips comments ('#') and blank lines
    next GFF_READ_LOOP if ($_ =~ m/^#.*$|^\s*$/);

    # parses each record line into a simple hash
    my %record = %{ readGFF ($_) };

    # skips records with non-unique matches
    next GFF_READ_LOOP if ($record{'score'} != 1);

    if ( scalar keys %HoH && $record{'start'} > max keys %HoH) {
        sort_and_count (\%HoH, \%total_count);
    }

    # grabs methylated sequence from gff file
    # assumes sequence is exactly $readsize bps long
    $record{'feature'} =~ m/([ACGTN]+$)/ or die Dumper \%record;
    my @methylated = split(//, $1);

    # grabs scaffold sequence from gff file
    # assumes sequence is in the 'attribute' field, separated from extraneous information by a '='
    my @unmethylated = split(//, (split "=", $record{'attribute'})[1]);

    # loops through each character in current sequences
    for ( my $i = $record{'start'}; $i < $record{'end'}; $i++ ) {
        
	# sets $j to 0 for forward strand coordinates
	my $j = $i - $record{'start'};

        # sets $k to $record{'end'} for reverse strand coordinates
	my $k = $record{'end'} - $j;

	# Since the input coordinates are absolute (ie. from the forward strand)
	# if a particular read maps to the reverse strand, the coordinates for any
	# 'c's hit need to be reversed
	my $coord = $i;
	$coord = $k if($record{'strand'} eq '-');

        # checks that we're looking at a left/1 sequence AND current character is a 'C'
	# this regex is adapted to the Solexa sequences we have, so it might change in the future
	# we're looking 2 characters ahead in the scaffold because it has 2 extra chars in the beginning
 	if ( $unmethylated[$j + 2] =~ m/[Cc]/
                 && $record{'feature'} =~ m/\/1:[ACGTN]+/ ) {

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

	    # grab some necessary information into the data structure
	    # we will print this information later
	    $HoH{$coord}[5] = $coord; # coord
	    $HoH{$coord}[6] = $record{'seqname'}; # chromosome
	    $HoH{$coord}[7] = $record{'strand'}; # strand
	}
	
	# checks that we're looking at a right/2 sequence AND current character is a 'C'
	# this regex is adapted to the Solexa sequences we have, so it might change in the future
	# we're looking 2 characters before in the scaffold because it has 2 extra chars in the beginning
	# AND we're looking in the reverse strand, so the context goes right-to-left
	elsif ( $unmethylated[$j + 2] =~ m/[Gg]/ && 
	      $record{'feature'} =~ m/\/2:[ACGTN]+/ ) {

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

	    # grab some necessary information into the data structure
	    # we will print this information later
	    $HoH{$coord}[5] = $coord; # coord
	    $HoH{$coord}[6] = $record{'seqname'}; # chromosome
	    $HoH{$coord}[7] = $record{'strand'}; # strand
	}
    }
}
close($GFF);
close(STDOUT);

open my $FREQ_OUT, '>', "$output.log" or die "Can't write to file: $output.log";

print $FREQ_OUT Dumper \%total_count;

close $FREQ_OUT;

# done

sub gff_sort {
    my $data_ref = shift;

    print STDERR "sorting data...\n";

    return sort {
        (split /\t/, $a)[3] <=> (split /\t/, $b)[3]
    } @{$data_ref};
}


sub sort_and_count {
    my ($single_count_ref, $total_count_ref) = @_;
    my %single_count = %{$single_count_ref};
    my %total_count = %{$total_count_ref};
    
    # loops through every initialized key in main hash
    # keys are sorted, so the output is going to be sorted by starting coordinate
    for my $i (sort {$a <=> $b} keys %single_count) { # Sorting and counting...  % done

        # we need to check that a given key/value pair was initialized
        # if it wasn't, initialize it to zero (to avoid division-by-zero, etc)
        for (0..4) {
            $single_count{$i}[$_] = 0 if !exists $single_count{$i}[$_];
        }

        $total_count{c_count} += $single_count{$i}[0];
        $total_count{t_count} += $single_count{$i}[1];
        $total_count{cg_count} += $single_count{$i}[2];
        $total_count{chg_count} += $single_count{$i}[3];
        $total_count{chh_count} += $single_count{$i}[4];
        %{$total_count_ref} = %total_count;       
        
        # calculates c/(c+t), which varies between 0 and 1
        # 0 is no methylation, 1 is total methylation
        if ($single_count{$i}[0] + $single_count{$i}[1] == 0) {
            $single_count{$i}[8] = 0;
        } else {
            $single_count{$i}[8] = $single_count{$i}[0]/($single_count{$i}[0]+$single_count{$i}[1]);
        }

        # finds the appropriate context to put in 'feature' field
        my $context = ".";
        if ($single_count{$i}[2]>$single_count{$i}[3] && $single_count{$i}[2]>$single_count{$i}[4]) {
            $context = "CG";
        } elsif ($single_count{$i}[4]>$single_count{$i}[3]) {
            $context = "CHH";
        } elsif ($single_count{$i}[3]>0) {
            $context = "CHG";
        }

        # checks that the context for any given 'c' is determinable
        # if not just output a '.' in its place
        if (($single_count{$i}[2] > 0 && $single_count{$i}[3] > 0) ||
       ($single_count{$i}[2] > 0 && $single_count{$i}[4] > 0) ||
       ($single_count{$i}[3] > 0 && $single_count{$i}[4] > 0)) {
            $context = ".";
        }

        # prints a single gff record
        print join("\t",
                   $single_count{$i}[6],
                   "dz_cm",
                   $context,
                   $single_count{$i}[5],
                   $single_count{$i}[5],
                   sprintf("%e", $single_count{$i}[8]), # score
                   $single_count{$i}[7],
                   ".",
                   join(";", "c=$single_count{$i}[0]", "t=$single_count{$i}[1]")
               ), "\n";

        delete $single_count{$i};
    }
    %{$single_count_ref} = ();
}


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
\t[--sort]\tPre-sort input file in-place before processing
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
