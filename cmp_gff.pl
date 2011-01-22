#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/sum/;

# Globals, passed as command line options
my $gff_file_1  = q{};
my $gff_file_2  = q{};
my $operation   = 'sub';
my $statistic   = 'Fisher::twotailed';
my $inverse_log = -1;
my $reverse     = 0;
my $min_meth    = 0;
my $min_diff    = 0;
my $min_sites   = 0;
my $skip_empty  = 0;
my $threshold   = 0;
my $concatenate = 0;
my $valid_gff   = 0;
my $ignore_feat = 0;
my $new_feature = 0;
my $debug       = 0;
my $output      = 0;
my $verbose     = 0;
my $quiet       = 0;
my $usage       = 0;

# Grabs and parses command line options
my $result = GetOptions (
    'gff-a|a=s'            => \$gff_file_1,
    'gff-b|b=s'            => \$gff_file_2,
    'operation|op|p:s'     => \$operation,
    'statistic|stat|s=s'   => \$statistic,
    'inverse-log|ilog|i:i' => \$inverse_log,
    'reverse-score|rev|r'  => \$reverse,
    'min-methylation|mm=f' => \$min_meth,
    'min-difference|md=f'  => \$min_diff,
    'min-sites|ms=i'       => \$min_sites,
    'skip-empty|k'         => \$skip_empty,
    'threshold|t'          => \$threshold,
    'concatenate|c'        => \$concatenate,
    'valid-gff|g'          => \$valid_gff,
    'ignore-features|if'   => \$ignore_feat,
    'new-feature|nf=s'     => \$new_feature,
    'debug|d'              => \$debug,
    'output|o=s'           => \$output,
    'verbose|v'            => sub { use diagnostics; },
    'quiet|q'              => sub { no warnings; },
    'help|h'               => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'             => sub { pod2usage ( -verbose => 2 ); }
);

my %statistics = (
    'CHI::phi'          => 1,
    'CHI::tscore'       => 1,
    'CHI::x2'           => 1,
    'Dice::dice'        => 1,
    'Dice::jaccard'     => 1,
    'Fisher::left'      => 1,
    'Fisher::right'     => 1,
    'Fisher::twotailed' => 1,
    'MI::ll'            => 1,
    'MI::pmi'           => 1,
    'MI::ps'            => 1,
    'MI::tmi'           => 1,
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $result and $gff_file_1 and $gff_file_2 and exists $statistics{$statistic};

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

# use the appropriate statistic measure based on user input
if (exists $statistics{$statistic}) {
    eval "use Text::NSP::Measures::2D::$statistic";
}

# opens gff files
open my $GFFA, '<', $gff_file_1 or croak "Can't read file: $gff_file_1";
open my $GFFB, '<', $gff_file_2 or croak "Can't read file: $gff_file_2";

# @window_buffers hold contiguous gff records, in the hash form produced by gff_read
#
# These buffers get filled while the current gff line being processed is contiguous 
# to the last one and get flushed (ie. processing resumes and the buffer is emptied) 
# when the current read is not contiguous to the last one
#
# The algorithm for determining contiguity is not very sophisticated: two windows are 
# adjacent or overlapping if the current line's starting coordinate is larger than the
# last processed line's start coordinate but smaller than, or equal to, the last 
# processed line's end coordinate
my (@window_buffer_a, @window_buffer_b) = ();

# processing stops when either we run out of 'a' file lines
# OR we run out of 'b' file lines (see below)
PROCESSING:
while (defined (my $line_a = <$GFFA>) and defined (my $line_b = <$GFFB>)) {

    # skip past comments on both files
    while ($line_a =~ m/^#.*$|^\s*$/) {
        $line_a = <$GFFA>;
    }
    while ($line_b =~ m/^#.*$|^\s*$/) {
        $line_b = <$GFFB>;
    }

    chomp $line_a;
    chomp $line_b;

    # read each line into a hash
    # see gff_read() sub definition for hash keys
    my %rec_a = %{&gff_read ($line_a)};
    my %rec_b = %{&gff_read ($line_b)};

    unless ($ignore_feat) {
        while ($rec_a{'feature'} ne $rec_b{'feature'}) { 
            if ($rec_a{'start'} > $rec_b{'start'}) {
                $line_a = <$GFFA>;
                chomp $line_a;
                %rec_a = %{&gff_read ($line_a)};
            } 
            elsif ($rec_a{'start'} < $rec_b{'start'}) {
                $line_b = <$GFFB>;
                chomp $line_b;
                %rec_b = %{&gff_read ($line_b)};
            }
        }
    }

    if ($skip_empty) {
        next PROCESSING
        if $rec_a{score} eq q{.}
        or $rec_b{score} eq q{.};
    }

    if ($min_sites) {
        my @a_sites = $rec_a{attribute} =~ m/[ct]=(\d+)/g;
        my @b_sites = $rec_b{attribute} =~ m/[ct]=(\d+)/g;

        next PROCESSING
        unless sum (@a_sites) >= $min_sites
        and    sum (@b_sites) >= $min_sites;
    }

    if ($min_meth) {
        next PROCESSING 
        unless $rec_a{score} >= $min_meth
        or     $rec_b{score} >= $min_meth;
    }

    if ($min_diff) {
        next PROCESSING
        unless abs($rec_a{score} - $rec_b{score}) >= $min_diff;
    }

    my $ngram = 0;
    if ($concatenate) {

        # filter out windows with no Cs
        next PROCESSING if
        ($rec_a{'attribute'} eq q{.} and $rec_b{'attribute'} eq q{.}) or
        ($rec_a{'attribute'} =~ m/c=0/ and $rec_b{'attribute'} =~ m/c=0/);

        # if our buffers are empty OR
        # if we think the current window is contiguous to the ones already stored
        # push the current windows into the buffers
        # and skip processing for now
        if (@window_buffer_a == 0) {
            push @window_buffer_a, $line_a;
            push @window_buffer_b, $line_b;
            next PROCESSING;
        }
        elsif ($rec_a{'start'} >  (split "\t", $window_buffer_a[-1])[3] and
               $rec_a{'start'} <= (split "\t", $window_buffer_a[-1])[4] + 1) {
            push @window_buffer_a, $line_a;
            push @window_buffer_b, $line_b;
            next PROCESSING;
        }
        else {
            # if the current window is NOT contiguous with the last buffer entry
            # we flush the buffer, start filling it up again, and overwrite the current
            # window with the concatenated window
            my $tmp_window_a_ref = gff_concatenate (\@window_buffer_a);
            my $tmp_window_b_ref = gff_concatenate (\@window_buffer_b);
            $ngram = gff_calculate_statistic ($tmp_window_a_ref, $tmp_window_b_ref);

            # deletes buffer contents, saves current windows to the buffer
            # and overwrites current windows and ngram score for printing
            (@window_buffer_a, @window_buffer_b) = ();
            push @window_buffer_a, $line_a;
            push @window_buffer_b, $line_b;

            # filter out records with statistic measures above threshold
            next PROCESSING if ($threshold and $ngram > $threshold);

            %rec_a = %{$tmp_window_a_ref};
            %rec_b = %{$tmp_window_b_ref};
        }

        # filter out windows with no Cs
        next PROCESSING if
        ($rec_a{'attribute'} eq q{.}   and $rec_b{'attribute'} eq q{.}) or
        ($rec_a{'attribute'} =~ m/c=0/ and $rec_b{'attribute'} =~ m/c=0/);
    }
    else {
        $ngram = gff_calculate_statistic (\%rec_a, \%rec_b);
    }

    my $score = 0;
    if ($operation eq 'sub') {
        $score = $rec_a{'score'} - $rec_b{'score'};
    }
    elsif ($operation eq 'div' and
           $rec_b{'score'} != 0) {
        $score = $rec_a{'score'} / $rec_b{'score'};
    }

    # calculates the inverse log (with base specified by user) of the ngram
    if ($ngram > 0 and $ngram != 1 and $inverse_log != -1) {
        if ($inverse_log == 0) {
            $ngram = 1 / log($ngram)
        }
        else {
            $ngram = 1 / (log($ngram) / log($inverse_log));
        }
    }

    # puts the score in the 'attribute' field and the statistic in the 'score' field
    if ($reverse) {
        my $tmp = $score;
        $score = $ngram;
        $ngram = $tmp;
        $statistic = $operation;
    }

    ### NOTE: this is temporary: it naively parses the input file names to try to
    # put something meaningful in the 'feature' field
    my $feature = $new_feature || "$rec_a{'feature'}:$operation:$rec_b{'feature'}";

    $statistic =~ tr/A-Z/a-z/;

    my $attribute = sprintf("$statistic=%g", $ngram);

    if ($debug) {
        $rec_a{'attribute'} =~ s/([ct]=)/a_$1/g;
        $rec_b{'attribute'} =~ s/([ct]=)/b_$1/g;
        $attribute .=  ';' . $rec_a{'attribute'} . ';' . $rec_b{'attribute'};
    }

    unless ($valid_gff) {
        my @attr_fields = (split /;/, $attribute);
        $attr_fields[$_] = (split /=/, $attr_fields[$_])[-1] for 0 .. $#attr_fields;
        $attribute = join "\t", @attr_fields;
    }

    # prints out the current window (or concatenated windows) as a gff record
    print join("\t",
               $rec_a{'seqname'},
               'cmp_gff.pl',
               $feature,
               $rec_a{'start'},
               $rec_a{'end'},
               sprintf("%g", $score),
               ".",
               ".",
               $attribute
           ), "\n";
}

close $GFFA;
close $GFFB;

exit 0;

=head1 Subroutines

=head2 gff_calculate_statistic()

    Takes two references to hashes output by gff_read().

    Returns an $ngram value determined by the choice of statistical measure

    Depends on the Text::NSP module by Ted Pedersen <http://search.cpan.org/dist/Text-NSP/>

=cut


sub gff_calculate_statistic {
    # unpack arguments
    my ($rec_a_ref, $rec_b_ref, $ignore_feature) = @_;
    my %rec_a = %{$rec_a_ref};
    my %rec_b = %{$rec_b_ref};
    $ignore_feature //= 1;

    # basic sanity test for matching coordinates, chromosome and context
    if ($rec_a{'start'} != $rec_b{'start'} or
        $rec_a{'end'} != $rec_b{'end'} or
        $rec_a{'seqname'} !~ m/$rec_b{seqname}/i or
        ($rec_a{'feature'} ne $rec_b{'feature'} && not $ignore_feature)) {
        print STDERR Dumper (\%rec_a, \%rec_b);
        croak (
            "Can't match windows in input files.
             Make sure these are sorted by starting
             coordinate and that each input file
             contains only one chromosome and one context.
             This will be modified in the future to allow
             multiple chromosomes and contexts per input file.",
        );
    }

    # checks for no coverage in the windows (ie. 'attribute' field eq ".")
    my $ngram = 0;
    if ( ($rec_a{'attribute'} =~ m/\./) or ($rec_b{'attribute'} =~ m/\./) ) {
        $ngram = 0;
    }
    else {
        # contigency table:
        #
        #  |line a|line b|
        # -|------|------|----
        # c| n11  | n12  | n1p
        # -|------|------|----
        # t| n21  | n22  |
        #  |------|------|----
        #  | np1  |      | npp 
        my ($n11) = $rec_a{'attribute'} =~ m/c=(\d+)/;
        my ($n21) = $rec_a{'attribute'} =~ m/t=(\d+)/;
        my ($n12) = $rec_b{'attribute'} =~ m/c=(\d+)/;
        my ($n22) = $rec_b{'attribute'} =~ m/t=(\d+)/;

        if ($n11 + $n12 + $n21 + $n22) {
            $ngram = calculateStatistic
            (
                n11 => $n11,
                n1p => ($n11 + $n12),
                np1 => ($n11 + $n21),
                npp => ($n11 + $n12 + $n21 + $n22)
            );
        }

        if ( my $error_code = getErrorCode() ) {
            print STDERR $error_code, ": ", getErrorMessage(), "\n";
        }
    }
    return $ngram;
}


=head2 gff_concatenate()

    Takes a reference to an array of contiguous windows in the form of gff_read() hashes

    Returns a similarly formatted hash with all contiguous windows merged and scores recalculated

=cut
sub gff_concatenate {
    # unpack argument
    my $contig_windows_ref = shift;
    my @contig_windows = @{$contig_windows_ref};

    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute);
    my ($c_count, $t_count) = (0, 0);
    # loops through every contiguous window
    for my $x (0..$#contig_windows) {

        my %rec = %{gff_read ($contig_windows[$x])};

        if ($x > 0) {
            # basic sanity test on whether current line is consistent with previous one
            if ($seqname  !~ m/$rec{'seqname'}/i
                or $source  ne $rec{'source'}
                or $feature ne $rec{'feature'}
                or $strand  ne $rec{'strand'}
                or $frame   ne $rec{'frame'}
                or $end + 1 <  $rec{'start'}) {
                croak ("Given records to be merged are inconsistent")
            }
        }
        # if ok, save current line
        else {
            $seqname = $rec{'seqname'};
            $source  = $rec{'source'};
            $feature = $rec{'feature'};
            $strand  = $rec{'strand'};
            $frame   = $rec{'frame'};
        }

        # get the lowest coordinate, naively assuming input windows are sorted
        # well, it is checked above, so we should never get an unsorted array
        $start = $rec{'start'} if ($x == 0);
        $end = $rec{'end'};
        $score += $rec{'score'} if $rec{'score'} ne q{.};

        # extract c and t counts from 'attribute' field
        # and update the total counts
        if ($rec{'attribute'} ne '.') {

            my ($non_over_c) = $rec{'attribute'} =~ m/c=(\d+)/;
            my ($non_over_t) = $rec{'attribute'} =~ m/t=(\d+)/;
            $c_count += $non_over_c;
            $t_count += $non_over_t;
        }
    }

    if ($c_count + $t_count == 0) {
        $score = 0;
    } 
    else {
        $score = $c_count/($c_count + $t_count);
    }

    # format 'score' to only output 3 decimal places
    $score = sprintf ("%g", $score);

    # format 'attribute' field for consistency here
    $attribute = "c=$c_count;t=$t_count";

    return {
        'seqname'   => $seqname,
        'source'    => $source,
        'feature'   => $feature,
        'start'     => $start,
        'end'       => $end,
        'score'     => $score,
        'strand'    => $strand,
        'frame'     => $frame,
        'attribute' => $attribute
    };
}


=head2 gff_read()

    Takes in a single GFF line string

    Returns a reference to an anonymous hash with keys equal to the GFF v3 spec

=cut
sub gff_read {
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute)
    = split /\t/, $_[0];

    return {
        'seqname'   => $seqname,
        'source'    => $source,
        'feature'   => $feature,
        'start'     => $start,
        'end'       => $end,
        'score'     => $score,
        'strand'    => $strand,
        'frame'     => $frame,
        'attribute' => $attribute,
    };
}

=head1 NAME

 cmp_gff.pl - Comparative enrichment analysis of two GFF files

=head1 SYNOPSIS

 cmp_gff.pl -a tissue_1.gff -b tissue_2.gff -p sub -s Fisher::twotailed -t 0.000001 -c -g -o a_b_comparison.gff

=head1 DESCRIPTION

=head1 OPTIONS

 cmp_gff.pl [OPTION]... -a [FILE A] -b [FILE B]

 -a,  --gff-a            first GFF alignment input file
 -b,  --gff-b            second GFF alignment input file
 -p,  --operation        arithmetic operation on scores from a and b ('sub' or 'div')
 -s,  --statistic        type of indendence/significance statistic to use
                             statistics options: (run 'perldoc Text::NSP' for more information)
                                 CHI::phi             Phi coefficient measure
                                 CHI::tscore          T-score measure of association
    		                 CHI::x2              Pearson's chi squared measure of association
    		                 Dice::dice           Dice coefficient
    		                 Dice::jaccard        Jaccard coefficient
    		                 Fisher::left         Left sided Fisher's exact test
    		                 Fisher::right        Right sided Fisher's exact test
    		                 Fisher::twotailed    Two-sided Fisher's exact test
    		                 MI::ll               Loglikelihood measure of association
    		                 MI::pmi              Pointwise Mutual Information
    		                 MI::ps               Poisson-Stirling measure of association
    		                 MI::tmi              True Mutual Information
 -i,  --inverse-log      compute inverse log of score in base specified by user
 -r,  --reverse-score    output scores in attributes field and statistics in scores field
 -t,  --threshold        maximum (optional) threshold for filtering out windows by p-value.
 -c,  --concatenate      concatenate adjacent windows after filtering stage
 -g,  --valid-gff        output format is valid GFF (attribute fields split by ';', not '\t'
 -d,  --debug            output individual c and t site counts
 -mm, --min-methylation  pre-filter by minimum methylation on either tissue
 -md, --min-difference   pre-filter by minimum difference on tissue comparison
 -ms, --min-sites        pre-filter by minimum total sites on each tissue
 -k,  --skip-empty       pre-filter by no coverage on either tissue
 -if, --ignore-features  don't match feature tracks
 -o,  --output           filename to write results to (defaults to STDOUT)
 -v,  --verbose          output perl's diagnostic and warning messages
 -q,  --quiet            supress perl's diagnostic and warning messages
 -h,  --help             print this information
 -m,  --manual           print the plain old documentation page

=head1 REVISION

 Version 0.0.2

 $Rev: 296 $:
 $Author: psilva $:
 $Date: 2010-04-05 16:40:53 -0700 (Mon, 05 Apr 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/cmp_gff.pl $:
 $Id: cmp_gff.pl 296 2010-04-05 23:40:53Z psilva $:

=head1 AUTHOR

 Pedro Silva <psilva@nature.berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
