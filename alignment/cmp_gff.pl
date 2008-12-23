#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use Carp;

# Globals, passed as command line options
my $gff_file_1 = '';
my $gff_file_2 = '';
my $operation = 'sub';
my $statistic = 'Fisher::right';
my $inverse_log = exp(1);
my $reverse = 0;
my $threshold = 0;
my $output = '-';
my $verbose = 0;
my $quiet = 0;
my $usage = 0;

# Initial check of command line parameters
if (@ARGV < 2) {usage();}
my @argv = @ARGV;

# Grabs and parses command line options
my $result = GetOptions ( 
    'gff-a|a=s' => \$gff_file_1,
    'gff-b|b=s' => \$gff_file_2,
    'operation|op|p:s' => \$operation,
    'statistic|stat|s=s' => \$statistic,
    "inverse-log|ilog|i:i" => \$inverse_log,
    'reverse-score|rev|r' => \$reverse,
    'output|o:s' => \$output,
    'threshold|t=f' => \$threshold,
    'verbose|v' => sub {enable diagnostics;use warnings;},
    'quiet|q' => sub {disable diagnostics;no warnings;},
    'usage|help|h' => \&usage
    );

# use the appropriate statistic measure based on user input
eval "use Text::NSP::Measures::2D::$statistic";

# redirects STDOUT to file if specified by user
if (!($output eq '-')) {
    open(STDOUT, '>', "$output") or die ("Can't redirect STDOUT to file: $output");
}

# opens gff files
open (my $GFFA, '<', $gff_file_1) or die ("Can't read file: $gff_file_1");
open (my $GFFB, '<', $gff_file_2) or die ("Can't read file: $gff_file_2");

# prints out header fields that contain gff v3 header, generating program, time, and field names
gff_print_header ($0, @argv);

my @contiguous_windows;
while (my $line_a = <$GFFA>) {

    chomp $line_a;
    next if ($line_a =~ m/^#.*$|^\s*$/);

    my $line_b = <$GFFB>;
    last unless defined $line_b;
    while ($line_b =~ m/^#.*$|^\s*$/) {
	$line_b = <$GFFB>;
    }

    $line_a =~ s/\r//g;
    $line_b =~ s/\r//g;

    my %rec_a = %{&gff_read ($line_a)};
    my %rec_b = %{&gff_read ($line_b)};
    
    ($line_a, my $context_a) = split(/;|_/, $rec_a{'feature'});
    ($line_b, my $context_b) = split(/;|_/, $rec_b{'feature'});

    if ($rec_a{'start'} != $rec_b{'start'} or
	$rec_a{'end'} != $rec_b{'end'} or
	$rec_a{'seqname'} ne $rec_b{'seqname'} or
	$context_a ne $context_b) {
	die ("Can't match windows in input files. Make sure these are sorted by starting coordinate\n
          and that each input file contains only one chromosome and one context. This will be\n
          modified in the future to allow multiple chromosomes and contexts per input file.");
    }
    
    my ($score, $ngram) = (0, 0);

    if ( ($rec_a{'attribute'} =~ m/\./) or
	 ($rec_b{'attribute'} =~ m/\./) ) {
	$ngram = -1;
	$score = 0;
    }
    else {
	my ($n11, $n21) = split(/;/, $rec_a{'attribute'});
	my ($n12, $n22) = split(/;/, $rec_b{'attribute'});

	($n11) = $n11 =~ m/(\d+)/;
	($n12) = $n12 =~ m/(\d+)/;
	($n21) = $n21 =~ m/(\d+)/;
	($n22) = $n22 =~ m/(\d+)/;
	
	$ngram = calculateStatistic (
	    n11 => $n11,
	    n1p => ($n11 + $n12),
	    np1 => ($n11 + $n21),
	    npp => ($n11 + $n12 + $n21 + $n22));

	if ( my $error_code = getErrorCode() ) {
	    print STDERR $error_code, ": ", getErrorMessage(), "\n";
	}
    }

    if ($threshold && $ngram >= $threshold) {


    }

    if ($operation eq 'sub') {$score = $rec_a{'score'} - $rec_b{'score'};}
    elsif ($operation eq 'div' && 
	   $rec_b{'score'} != 0) {$score = $rec_a{'score'} / $rec_b{'score'};}

    if ($ngram > 0 && $ngram != 1) {
	if ($inverse_log) {$ngram = 1 / (log($ngram) / log($inverse_log));}
	else {$ngram = 1 / log($ngram);}
    }

    if($reverse) {
	my $tmp = $score;
	$score = $ngram;
	$ngram = $tmp;
	$statistic = $operation;
    }

    print join("\t",
	       $rec_a{'seqname'},
	       "cmp",
	       "${line_a}${operation}${line_b};$context_a",
	       $rec_a{'start'},
	       $rec_a{'end'},
	       sprintf("%.3f", $score),
	       ".",
	       ".",
	       sprintf("$statistic=%.3f", $ngram)), "\n";
    }
    
close ($GFFA);
close ($GFFB);

exit 0;


sub gff_concatenate {
    my $contig_windows_ref = @_;
    my @contig_windows = @{$contig_windows_ref};

    my ($start, $end, $score, $c_count, $t_count, $seqname, $source, $feature, $strand, $frame);
    foreach my $line (@contig_windows) {
	my %rec = %{gff_read ($line)};
	$seqname = $rec{'seqname'};
	$source = $rec{'source'};
	$feature = $rec{'feature'};
	$strand = $rec{'strand'};
	$frame = $rec{'frame'};

	$start = $rec{'start'} if ($. == 0);
	$end = $rec{'end'} if ($. == scalar @contig_windows - 1);
	$score += $rec{'score'};

	my ($tmp1, $tmp2) = split(/;/, $rec{'attribute'});
	($c_count) += $tmp1 =~ m/(\d+)/;
	($t_count) += $tmp2 =~ m/(\d+)/;
    }
    $score /= scalar @contig_windows;

    return (
	'seqname' => $seqname,
	'source' => $source,
	'feature' => $feature,
	'start' => $start,
	'end' => $end,
	'score' => $score,
	'strand' => $strand,
	'frame' => $frame,
	'attribute' => "c=$c_count;t=$t_count"
	);
}



# readGFF reads in a single GFF line, and outputs a hash reference
sub gff_read {
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $_[0]);
    my %rec = (
	'seqname'=>$seqname,
	'source'=>$source,
	'feature'=>$feature,
	'start'=>$start,
	'end'=>$end,
	'score'=>$score,
	'strand'=>$strand,
	'frame'=>$frame,
	'attribute'=>$attribute
	);
    return \%rec;
}

# prints out a commented line with header fields
sub gff_print_header {
    print "##gff-version 3\n";
    print join(" ",
	       "#",
	       @_,
	       "\n");
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime (time);
    printf "# %4d-%02d-%02d %02d:%02d:%02d\n", $year+1900, $mon+1, $mday, $hour, $min, $sec;
    print join("\t",
	       "# SEQNAME",
	       "SOURCE",
	       "FEATURE",
	       "START",
	       "END",
	       "SCORE",
	       "STRAND",
	       "FRAME",
	       "ATTRIBUTES",
	), "\n";
}


# prints out usage information
sub usage {
    print STDERR <<EOF;
countMethylation.pl \<REQUIRED\> [OPTIONS]
    <--gff-a          -a>    First GFF alignment input file
    <--gff-b          -b>    Second GFF alignment input file
    [--operation      -p]    Arithmetic operation on scores from a and b (\'sub\' or \'div\')
    [--statistic      -s]    Type of indendence/significance statistic to use
                      Statistics options: (run \'perldoc Text::NSP\' for more information)
                      CHI::phi                 Phi coefficient measure 
		      CHI::tscore              T-score measure of association 
		      CHI::x2                  Pearson\'s chi squared measure of association 
		      Dice::dice               Dice coefficient 
		      Dice::jaccard            Jaccard coefficient
		      Fisher::left             Left sided Fisher\'s exact test
		      Fisher::right            Right sided Fisher\'s exact test
		      Fisher::twotailed        Two-sided Fisher\'s exact test
		      MI::ll                   Loglikelihood measure of association 
		      MI::pmi                  Pointwise Mutual Information
		      MI::ps                   Poisson-Stirling measure of association
		      MI::tmi                  True Mutual Information
    [--inverse-log    -i]    Step interval of sliding window in bp
    [--reverse-score  -r]    Output scores in attributes field and statistics in scores field 
    [--threshold      -t]    Minimum threshold for filtering out windows
    [--output         -o]    Filename to write results to (default is STDOUT)
    [--verbose        -v]    Output perl\'s diagnostic and warning messages
    [--quiet          -q]    Supress perl\'s diagnostic and warning messages
    [--usage          -h]    Print this information
EOF
    exit 0;
}
