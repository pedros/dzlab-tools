#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use Data::Dumper;
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

my (@window_buffer_a, @window_buffer_b) = ();
while (my $line_a = <$GFFA>) {

    next if ($line_a =~ m/^#.*$|^\s*$/);
    my $line_b = <$GFFB>;
    last unless defined $line_b;
    while ($line_b =~ m/^#.*$|^\s*$/) {$line_b = <$GFFB>}
    chomp $line_a;
    chomp $line_b;

    $line_a =~ s/\r//g;
    $line_b =~ s/\r//g;

    my %rec_a = %{&gff_read ($line_a)};
    my %rec_b = %{&gff_read ($line_b)};

    my $ngram = gff_calculate_statistic (\%rec_a, \%rec_b);
    my $score = 0 if $ngram == -200;

     if ($threshold) {

	 next if ($ngram < $threshold);
	
	 if (@window_buffer_a == 0 or
             (@window_buffer_a > 0 and 
              $rec_a{'start'} > $window_buffer_a[$#window_buffer_a]{'start'} and 
              $rec_a{'start'} <= $window_buffer_b[$#window_buffer_b]{'end'} + 1)) {
		 push @window_buffer_a, %rec_a;
		 push @window_buffer_b, %rec_b;
                 next;
         }
         else {
             # flush the buffer and start filling it up again
             my $tmp_window_a = gff_concatenate (@window_buffer_a);
             my $tmp_window_b = gff_concatenate (@window_buffer_b);
             my $tmp_ngram = gff_calculate_statistic ($tmp_window_a, $tmp_window_b);
             (@window_buffer_a, @window_buffer_b) = ();
             push @window_buffer_a, %rec_a;
             push @window_buffer_b, %rec_b;
             %rec_a = %{$tmp_window_a};
             %rec_b = %{$tmp_window_b};
             $ngram = $tmp_ngram;
         }
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

    my $tmp = substr($gff_file_1, 5, 3) . $operation . substr($gff_file_2, 5, 3);

    print join("\t",
	       $rec_a{'seqname'},
	       "cmp",
	       $tmp,
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


sub gff_calculate_statistic {
    my ($rec_a_ref, $rec_b_ref) = @_;
    my %rec_a = %{$rec_a_ref};
    my %rec_b = %{$rec_b_ref};

    my ($line_a, $context_a) = split(/;|_/, $rec_a{'feature'});
    my ($line_b, $context_b) = split(/;|_/, $rec_b{'feature'});

    if ($rec_a{'start'} != $rec_b{'start'} or
	$rec_a{'end'} != $rec_b{'end'} or
	$rec_a{'seqname'} ne $rec_b{'seqname'} or
	$context_a ne $context_b) {
	die ("Can't match windows in input files. Make sure these are sorted by starting coordinate\n
          and that each input file contains only one chromosome and one context. This will be\n
          modified in the future to allow multiple chromosomes and contexts per input file.");
    }
    
    my $ngram = 0;
    if ( ($rec_a{'attribute'} =~ m/\./) or
	 ($rec_b{'attribute'} =~ m/\./) ) {
	$ngram = -200;
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
    return $ngram;
}


sub gff_concatenate {
    my $contig_windows_ref = $_[0];
    my @contig_windows = @{$contig_windows_ref};

    my ($start, $end, $score, $seqname, $source, $feature, $strand, $frame, $last_end);
    my ($c_count, $t_count) = (0, 0);
    for my $x (0..$#contig_windows) {

# 	my $line = $contig_windows[$x];
# 	chomp $line;
# 	my %rec = %{gff_read ($line)};	

	my %rec = $contig_windows[$x];

	if  ($x > 0) {
	    if ($seqname ne $rec{'seqname'}
		or $source ne $rec{'source'}
		or $feature ne $rec{'feature'}
		or $strand ne $rec{'strand'}
		or $frame ne $rec{'frame'}
		or ($end > $rec{'start'} and $end < $rec{'end'}))
	    {
		die ("Given records to be merged are inconsistent")
	    }
	}
	else {
	    $seqname = $rec{'seqname'};
	    $source = $rec{'source'};
	    $feature = $rec{'feature'};
	    $strand = $rec{'strand'};
	    $frame = $rec{'frame'};
	}

	$start = $rec{'start'} if ($x == 0);
	$end = $rec{'end'};
	$score += $rec{'score'};

	my ($tmp1, $tmp2) = split(/;/, $rec{'attribute'});
	($tmp1) = $tmp1 =~ m/(\d+)/;
	($tmp2) = $tmp2 =~ m/(\d+)/;
	$c_count += $tmp1;
	$t_count += $tmp2;
    }

    $score = sprintf ("%.3f", ($c_count/($c_count+$t_count)));

    my %read_hash = (
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
    return \%read_hash;
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
