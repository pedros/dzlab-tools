#!/usr/bin/perl

=head1 NAME
chipotle.pl - find significant features in tiling array ChIP-chip data

=head1 SYNOPSIS

  % chipotle.pl --infile FILE.dat --model gauss --window 1000 --step 150 --correction BH --alpha 0.05 --dolog2 --transform 0.0
  

=head1 DESCRIPTION
Perl port of Michael Buck's ChIPotle algorithm (PMID: 16277752). It currently supports only the Gaussian background model method for significance estimation, while Buck's implementation can also use permutation analysis. Multiple testing correction is implemented using the Benjamini/Hochberg (1995) False Discovery Rate method.

=head1 COMMAND-LINE OPTIONS

Command-line options cannot be abbreviated to single-letter options.
   --infile     <filename>    File containing feature coordinates and values
   --model      <model>       Significance estimation mode (default 'gauss')
   --window     <integer>     Sliding window size in base pairs (default 990)
   --step       <integer>     Step size in base pairs (default 198)
   --transform  <float>       Non-zero array mean correction factor (default 0)
   --correction <float>       Multiple testing correction method (default 'BH')
   --alpha      <float>       Alpha for multiple testing correction (default 0.05)
   --dolog2                   Log2 transform array values (default is false)
   --usage                    Print usage information and quit


=head1 INPUT FILE FORMAT
In the original implementation, microarray data is pasted into an Excel worksheet. Columns are then selected via a macro. The present implementation expects an tab-delimited input file with specific columns. Data should be transformed into this file format using other scripts that work directly on hybridization files such as .GPR and .PAIR files or by querying a microarray database.
The columns expected, in order in a tab-delimited file, are:
ElementID: The name of a given microarray probe
Reference sequence: This is the reference sequence for a given array feature. Examples are genbank accessions, chromosome IDs, and BAC numbers.
Start position: The genomic start location for a given microarray probe
End position: The genomic end location for a given microarray probe
Log2(ratio): If your data is not log-transformed, this version of Chipotle can actually do that for you if you pass the --dolog2 flag on the command line
=cut


=head1 AUTHOR
Matthew Vaughn, vaughn@cshl.edu
Copyright (c) 2006 Cold Spring Harbor Laboratory

=head1 SUBROUTINES

=cut

#use strict;
use warnings;
use POSIX qw(floor);
use Statistics::Descriptive::Discrete;
use Getopt::Long;
# Globals, passed as command line options
my $infile   = "experiment.txt"; # input file name
my $model = 'gauss';             # significance estimation model
my $windowSize = 300;            # window size
my $stepSize = 38;         # how far to slide window in each iteration
my $doLogTransform = 0;    # log2 transform ratios in file
my $alpha = 0.05;          # alpha for multiple testing correction
my $method = "BH";         # method for MTC
my $transform = 0; # shift all values by this amount. Used to correct for non-zero array mean
my $usage = 0;     # print usage and exit
my $result = GetOptions (	"infile=s" => \$infile,
                                "model=s" => \$model,
                                "window=i"   => \$windowSize,
                                "step=i" => \$stepSize,
                                "correction=s" => \$method,
                                "alpha=f" => \$alpha,
                                "dolog2"  => \$doLogTransform,
                                "transform=f" => \$transform,
                                "usage|help" => \$usage);
if ($usage) {
    print STDOUT "chipotle.pl --infile filename [experiment.txt] --model string(gauss,permutation,none) --window int [1000] --step int [250] --correction string(None,BON,BH) [BH] --alpha float [0.05] --dolog2 flag [F] --usage", "\n";
    exit 1;
}
# These defined the column IDs for the input file
use constant kElementCol => 3;
use constant kReferenceCol => 0;
use constant kStartCol => 3;
use constant kEndCol => 4;
use constant kLogRatioCol => 5;
# Columns for the window-based output file
use constant kOutWindowIDCol => 0;
use constant kOutReferenceCol => 1;
use constant kOutStartCol => 2;
use constant kOutEndCol => 3;
use constant kOutValueCol => 4;
use constant kOutPvalueCol => 5;
use constant kOutSigCallCol => 6;
# States for initial data-reading Finite State Machine
use constant kGenericFailure => -1;
use constant kReadLine => 0;
use constant kWithinWindow => 1;
use constant kWrongReferenceSequence => 2;
use constant kPastWindow => 3;
use constant kBeforeWindow => 4;
use constant kNewWindow => 5;
use constant kWithinFeature => 6;

# Opening a file to capture standard output:  infile_alpha.STDOUT
my $oOUT = $infile;
$oOUT =~ s/\..*/\_win$windowSize\_step$stepSize\_alpha$alpha\.STDOUT/g;
open(SOUT, ">$oOUT");
print SOUT "Processing data in $infile\n";
print SOUT "Using $windowSize bp window and $stepSize bp steps\n";
print STDERR ".......Sending STDOUT for ChIPOTle to $oOUT\n";
# Pass through data once and grab the background features
# This is farmed out to a subroutine to keep the main logic clean
my $referenceSequenceBegin;
my $backgroundStdev = CalculateBackgroundStdev($infile);
print SOUT "Background Stdev = ", sprintf("%.3f", $backgroundStdev), "\n";

# Pass through data again and calculate signal and p-value for each window
# These arrays are globals used to store results from the SlidingWindow sub
my (@windowDataRows, @pValueCollection);
my $numWindows = RunSlidingWindow($infile);

# Now perform multiple testing correction
my @windowDataRowsTested;
my $multipleTestingSignificanceThreshold = FindMultipleTestingCorrectionThreshold(\@pValueCollection, $alpha, $method);
print SOUT "For alpha = $alpha, MTCST: $multipleTestingSignificanceThreshold\n";
my $numSigWindows = EvaluateSlidingWindowResults();
print SOUT "After multiple testing correction, there were $numSigWindows significant windows\n";

# Now, slide through results and find collapse significant overlapping windows into significant peaks or troughs.
# Results are collected in yet another array of arrayrefs
my (@significantPeaks, @significantValleys);
my ($countPeaks, $countValleys) = CollapseWindowsIntoFeatures();

print SOUT "This translates into $countValleys significant negative peaks and $countPeaks positive peaks\n";
print STDERR "There were $countPeaks peaks\n";
print STDERR "Finished!\n";
close SOUT;
# END MAIN CODE

sub CollapseWindowsIntoFeatures {

=head2 CollapseWindowsIntoFeatures

Reads in intermediate features file and collapses significant, overlapping windows into peaks. Returns array containing count of peaks and valleys identified.

=cut
    # Note: This does not handle changes in reference sequence context.
    print SOUT "Collapsing significant, overlapping windows into peaks and valleys.\n";
    my ($flagGAP, $countPeaks, $countValleys, $collectingFeature) = (0,0,0,0);
    my ($featureLeft, $featureRight, $featureID, $featureRefseq) = (0,0,0,undef);
    my $featureValStat = Statistics::Descriptive::Discrete->new();
    my $featurePvalStat = Statistics::Descriptive::Discrete->new();
    my $featureMaxWin = 0;
    # Opening peak output file
    my $out = $infile;
    $out =~ s/\..*/\_win$windowSize\_step$stepSize\_alpha$alpha\.peaks/g;
    open (FEATURES, ">$out") or die "Can't open output file for peaks: $out\n";
    foreach my $window (@windowDataRowsTested) {
        my @cols = @{ $window };
        if ($cols[kOutSigCallCol] == 1) {
            # Were we already in an established feature?
            if ($collectingFeature == 1) {
				# Yes, so extend its right coordinate to end of current window
				# Mapping Peaks from the Center of Windows
                $newRight = int( ($cols[kOutEndCol] + $cols[kOutStartCol]) /2 );
                $distChanged = $newRight - $featureRight;
                if ( $distChanged > $stepSize || $featureRefseq ne $cols[kOutReferenceCol]) {
                    $flagGAP = 1;
                } elsif ($distChanged == $stepSize) {
                    $featureRight = $newRight;
                    # and collect its value and p-value
                    $featureValStat->add_data($cols[kOutValueCol]);
                    $featurePvalStat->add_data($cols[kOutPvalueCol]);
                    $collectingFeature = 1;
                    # Updating maximum window value
                    if ($cols[kOutValueCol] > $featureMaxWin) {
                        $featureMaxWin = $cols[kOutValueCol];
                    }
                    $flagGAP = 0;
                }
                if ($flagGAP) {
                    my $value = sprintf("%.3f", $featureValStat->mean());
                    my $pvalue = sprintf("%.3e", $featurePvalStat->mean());
                    my @temp = ($featureID, $featureRefseq, $featureLeft, $featureRight, $value, $pvalue, $featureMaxWin);
                    print FEATURES join("\t", @temp), "\n";
                    if ($value > 0) {
                        $countPeaks++;
                    } elsif ($value < 0) {
                        $countValleys++;
                    }
                    $featureLeft = $newRight;
                    $featureRight = $newRight;
                    $featureRefseq = $cols[kOutReferenceCol];
                    # Keeping track of maximum window value
                    $featureMaxWin = $cols[kOutValueCol];
                    $featureID++;
                    $featureValStat = Statistics::Descriptive::Discrete->new();
                    $featurePvalStat = Statistics::Descriptive::Discrete->new();
                    $featureValStat->add_data($cols[kOutValueCol]);
                    $featurePvalStat->add_data($cols[kOutPvalueCol]);
                    $collectingFeature = 1;
                    $flagGAP = 0;
                }
            } else {
				# Establish a new feature and start adding to it
				# Mapping Peaks from the Center of Windows
                $featureLeft = int( ($cols[kOutEndCol] + $cols[kOutStartCol]) /2 );
                $featureRight = int( ($cols[kOutEndCol] + $cols[kOutStartCol]) /2 );
                $featureRefseq = $cols[kOutReferenceCol];
				# Keeping track of maximum window value
                $featureMaxWin = $cols[kOutValueCol];
                $featureID++;
                $featureValStat = Statistics::Descriptive::Discrete->new();
                $featurePvalStat = Statistics::Descriptive::Discrete->new();
                $featureValStat->add_data($cols[kOutValueCol]);
                $featurePvalStat->add_data($cols[kOutPvalueCol]);
                $collectingFeature = 1;
            }
        } else {
            if ($collectingFeature == 1) {
				# If the previous data row was in a feature and this row is not signficant, then we're not collecting windows into that feature
				# output that feature
                my $value = sprintf("%.3f", $featureValStat->mean());
                my $pvalue = sprintf("%.3e", $featurePvalStat->mean());
                if ($featureLeft == $featureRight) {
                    $extend = int($stepSize/2);
                    $featureLeft = $featureLeft - $extend;
                    $featureRight = $featureRight + $extend;
                }
                my @temp = ($featureID, $featureRefseq, $featureLeft, $featureRight, $value, $pvalue, $featureMaxWin);
                print FEATURES join("\t", @temp), "\n";
                $collectingFeature = 0;
                if ($value > 0) {
                    $countPeaks++;
                } elsif ($value < 0) {
                    $countValleys++;
                }
            }
        }
    }
    close FEATURES;
    return ($countPeaks, $countValleys);
}

sub EvaluateSlidingWindowResults {

=head2 EvaluateSlidingWindowResults

    This subroutine scans through the list of windows and compares their p-value to the multiple-testing-correction derived significance threshold. It returns the integer number of significanct features.

=cut

    print STDERR "Evaluating per-window results for significance.\n";
    my $significantFeatures = 0;
    #################################################################################################################
    # PRINTS VALUES FOR EACH WINDOW
    #	$out = $infile;
    #	$out =~ s/\..*/\_windows\.tsv/g;
    #	open (WINDOWS, ">$out") or die;
    foreach my $wref (@windowDataRows) {
        my @row = @{$wref};
        my $significant;
        if ($row[kOutPvalueCol] <= $multipleTestingSignificanceThreshold) {
            $significant = 1;
            $significantFeatures++;
        } else {
            $significant = 0;
        }
        push (@row, $significant);
        #		print WINDOWS join("\t", @row), "\n";
        push(@windowDataRowsTested, \@row);
    }
    #	close WINDOWS;
    #################################################################################################################
    return $significantFeatures;
}

sub RunSlidingWindow {

=head2 RunSlidingWindow

    The base code ported over from the ChIPotle algorithm. Evaluates a sliding window of genomic tiling ChIP data against a Gaussian error function to find significant peaks. Returns the integer ID of the last window examined, which is a proxy for the number of windows examined.

=cut

    print STDERR "\nRunning primary sliding window.\n";
    if ($doLogTransform) {
        print SOUT "Caution: Using on-the-fly log2 transformation\n";
    }
    my $filename = shift;
    open (IN, $filename) or die "Couldn't access $filename\n";
    my @rows = <IN>;
    # Filter out comments and such using grep
    @rows = grep { !/^\#/ } @rows;
    my @line = split(/\t/,$rows[0]);
    # Attributes of the current sliding window
    # References the first unique region
    my $windowLeft = $referenceSequenceBegin{$line[1]};
    my $windowRight = $windowLeft + $windowSize;
    my $windowCurrentReferenceSeq;
    # currentRow is the index within the array of data lines
    my ($currentRow, $windowCountFeatures, $windowID) = (0, 0, 1);	
    my $windowStatObj = Statistics::Descriptive::Discrete->new();
    my $finalRow = $#rows;
    print SOUT "There are ". ($finalRow + 1) . " rows of data in $filename\n";
    while ($currentRow <= $finalRow) {
        my $rowData = $rows[$currentRow];
        chomp($rowData);
        my @cols = split(/\s+/, $rowData);
        my $featurePos = floor(($cols[kEndCol] + $cols[kStartCol]) / 2);
        # Initialize current reference sequence
        if (! defined($windowCurrentReferenceSeq)) {
            $windowCurrentReferenceSeq = $cols[kReferenceCol];
        }
        # Establish status of current data point with respect to current sliding window
        # This is essentially a Finite-State Machine, albeit a poorly-implemented one
        my $currentState = kReadLine;
        if (($featurePos >= $windowLeft) and ($featurePos <= $windowRight) and ($cols[kReferenceCol] eq $windowCurrentReferenceSeq)) {
            $currentState = kWithinWindow;
        } elsif ($cols[kReferenceCol] ne $windowCurrentReferenceSeq) {
            $currentState = kWrongReferenceSequence;
        } elsif ($featurePos > $windowRight) {
            $currentState = kPastWindow;
        } elsif ($featurePos < $windowLeft) {
            $currentState = kBeforeWindow;
        } else {
            $currentState = kGenericFailure;
        }
		
        # Extract the value for the current data point, transforming if necessary
        my $value = $cols[kLogRatioCol];
        if ($doLogTransform) {
            # Protect from zero values
            if ($value <= 0) {
                $value = 0;
            } else {
                $value = log($value) / log(2);
            }
        }
        $value = $value - $transform;
        if ($currentState == kWithinWindow) {
            # Cache current reference sequence. Not actually needed, but formally we should do this while inside this loop
            $windowCurrentReferenceSeq = $cols[kReferenceCol];
            # As long as feature is within current window, keep extracting values
            $windowStatObj->add_data($value);
            $windowCountFeatures++;
            $currentRow++;
        } else {
            # Only compute window values if:
            # 1. We've run downstream of the current window
            # 2. There are one or more features in the array of values. 
            # The 2nd condition accounts for cases where a window cannot be filled. Primarily, this will be due to gaps in the genomic tiling array coverage.
            if (($currentState == kPastWindow) and ($windowCountFeatures > 0)) {
				# Calculate average ratio
                my $windowAverageLogRatio = $windowStatObj->mean();
				# Calculate uncorrected pvalue
                my $windowPvalue = 1;
                $windowPvalue = ErrorFunction( $windowAverageLogRatio / ($backgroundStdev / sqrt($windowCountFeatures)) ) / 2;
				# Output to temporary file
                my @temp = ($windowID, $cols[kReferenceCol], $windowLeft, $windowRight, $windowAverageLogRatio, $windowPvalue);
				# Store references to the window data
                push(@windowDataRows, \@temp);
				# Store all p-values so I can do multiple testing correction on the entire data set later
                push(@pValueCollection, $windowPvalue);			
            }
            # Slide the window N base pairs from start of previous one if and only if we have entered this loop due to passing the end of the previous window
            if ($currentState == kPastWindow) {
				# Set up the sliding window's new coordinates
                $windowLeft = $windowLeft + $stepSize;
				# Backup until feature is within new window. This is not robust for dealing with array gaps or changes in reference sequence
                if ($featurePos > $windowLeft) {
                    my $tempPos = floor(($cols[kEndCol] + $cols[kStartCol]) / 2);
                    while (($featurePos > $windowLeft) and ($currentRow > 0)) {
                        $currentRow--;
                        $rowData = $rows[$currentRow];
                        unless ($rowData =~ /^[#!\s]/) {
                            chomp($rowData);
                            my @cols = split(/\s+/, $rowData);
                            $featurePos = floor(($cols[kEndCol] + $cols[kStartCol]) / 2);
                        }
                    }
                    $currentRow++;
                }
                $currentState = kNewWindow;
            }
            # Start the window over at the beginning of a new reference sequence if we've changed refseqs
            if ($currentState == kWrongReferenceSequence) {
				# Start at beginning of new reference sequence
                $windowLeft = $referenceSequenceBegin{$cols[kReferenceCol]};
                $currentState = kNewWindow;
            }
            if ($currentState == kNewWindow) {
                $windowRight = $windowLeft + $windowSize;
                $windowID++;
                $windowCurrentReferenceSeq = $cols[kReferenceCol];
				# Re-initialize Statistics::Descriptive for holding window data
                $windowStatObj = Statistics::Descriptive::Discrete->new();
                $windowCountFeatures = 0;
            }
            #			if ($currentState == kBeforeWindow) {$currentRow++;}
        }
    }
    return $windowID;
}

sub CalculateBackgroundStdev {

=head2 CalculateBackgroundStdev

This routine calculates the standard deviation of the negative-valued (unbound) portion of the ChIP array. In the original implementation, the entire array is loaded into memory to find the background distribution. That's not viable for a high-volume implemenation, so instead, I pass once through the entire data set and harvest all the background values.

=cut

    	print STDERR "\nCalculating background standard deviation\n";
	if ($doLogTransform) {
		print SOUT "Caution: Using on-the-fly log2 transformation\n";
	}
	my $filename = shift;
	# Keep count of background elements
	# Just a temporary value to hold the deviations
	# The standard deviation for background features across the entire data set
	my ($countBackgroundElements, $tempBackgroundStdev);
	my ($sumDeviation, $sumAllValues, $countAllElements, $regionOLD) = (0,0,0,'NA');
	open (IN, $filename) or die "Couldn't access $filename\n";
	while (my $row = <IN>) {
		# Ignore comments and whitespace lines
		unless ($row =~ /^[#!\s]/) {
			chomp($row);
			my @cols = split(/\s+/, $row);
			my $value = $cols[kLogRatioCol];
			# Captures the start (lowest) coordinate for each unique region (assume data is sorted)
			my $region = $cols[kReferenceCol];
			my $start = $cols[kStartCol];
			unless($regionOLD eq $region){$referenceSequenceBegin{$region} = 99999999999999999;}
			if($referenceSequenceBegin{$region} > $start){
				$referenceSequenceBegin{$region} = $start;
				print SOUT "\nSet $region start to $referenceSequenceBegin{$region}\n";
			}
			$regionOLD = $region;
			if ($doLogTransform) {
				# Protect from zero values
				if ($value <= 0) {$value = 0;}
				else {$value = log($value) / log(2);}
			}
			$value = $value - $transform;
			if ($value < 0) {
################################################################################################
# Subtracts each value from 0, no estimate of the mean or median
################################################################################################
				$sumDeviation = $sumDeviation + ((0 - $value) ** 2);
				$countBackgroundElements++;
			}
			$sumAllValues = $sumAllValues + $value;
			$countAllElements++;
		}
	}close IN;
	my $datasetMeanValue = sprintf("%.3f", ($sumAllValues / $countAllElements));
	print SOUT "Data set mean: $datasetMeanValue\n";
	if (abs(0 - $datasetMeanValue) > 0.5) {
		print STDERR "Warning: this data set's mean value is really skewed away from a mean of zero!\nRe-run chipotle using '--transform $datasetMeanValue'\n";
		exit 1;
	}
	$tempBackgroundStdev = sqrt($sumDeviation / ($countBackgroundElements - 1));
	return $tempBackgroundStdev;
}

sub FindMultipleTestingCorrectionThreshold {

=head2 FindMultipleTestingCorrectionThreshold

usage: FindMultipleTestingCorrectionThreshold(\@array_of_pvalues, 0.05, 'BH')
This subroutine takes an array of p-values, an alpha, and a MTC method name and returns a corrected significance threhold. Currently, only the Benjamini and Hochberg (1995) FDR method is supported. The subroutine returns the threshold as a floating point value.

=cut
	print STDERR "Finding significance threshold for multiple testing correction.\n";
	my $pvalArrayReference = shift;
	my $alpha = shift;
	my $method = shift;
	my $result = $alpha;
	if ($method eq 'BH') {
		$result = _CutoffBenjaminiAndHochberg1995($pvalArrayReference, $alpha);
	}elsif ($method eq 'BON') {
		$result = _CutoffBonferroni($pvalArrayReference, $alpha);
	}
	return $result;
}

sub _CutoffBonferroni {
	print SOUT "Method: Bonferroni Correction.\n";
	my $ref = shift;
	my $alpha = shift;
	my @unsortedPvalues = @{ $ref };
	# @sortedPvalues is a vector of p-values, sorted in ascending order
	my @sortedPvalues = sort { $a <=> $b } @unsortedPvalues;
	my $countPvalues = scalar @sortedPvalues;
	my $jAlpha = ($alpha / $countPvalues);
	my $diffPJ;
	my $thresholdIndex = 0;
	for my $d (0..($countPvalues - 1)) {
		$diffPJ = $sortedPvalues[$d] - $jAlpha;
		if ($diffPJ < 0) {$thresholdIndex = $d;}
	}
	my $threshold = $sortedPvalues[$thresholdIndex];
	return $threshold;
}

sub _CutoffBenjaminiAndHochberg1995 {

=head2 _CutoffBenjaminiAndHochberg1995

Private method: implementation of False Discovery Rate of Benjamini and Hochberg, 1995

=cut
	# http://www.unt.edu/benchmarks/archives/2002/april02/rss.htm
	print SOUT "Method: False Discovery Rate of Benjamini and Hochberg, 1995.\n";
	my $ref = shift;
	my $alpha = shift;
	my @unsortedPvalues = @{ $ref };
	# @sortedPvalues is a vector of p-values, sorted in ascending order
	my @sortedPvalues = sort { $a <=> $b } @unsortedPvalues;
	my $countPvalues = scalar @sortedPvalues;
	my @jAlpha;
	for my $j (1..$countPvalues) {
		my $bit = $j * ($alpha / $countPvalues);
		$jAlpha[$j-1] = $bit;
	}
	# @jAlpha is a series containing the same number of values as the number of P-values in the initial vector
	my @diffPJ;
	my $thresholdIndex = 0;
	for my $d (0..($countPvalues - 1)) {
		$diffPJ[$d] = $sortedPvalues[$d] - $jAlpha[$d];
		if ($diffPJ[$d] < 0) {$thresholdIndex = $d}
	}
	my $threshold = $sortedPvalues[$thresholdIndex];
	return $threshold;
}

sub GetMax {

=head2 GetMax

This is a subroutine name from the original ChIPotle source code. It returns scalar maximum value for array. It assumes all values in array are numeric.

=cut
	my $arrayref = shift or die;
	my @unsorted = @{ $arrayref };
	my @sorted = sort { $a <=> $b } @unsorted;
	my $max_index = $#sorted;
	return $sorted[$max_index];
}

sub GetMin {

=head2 GetMin

This is a subroutine name from the original ChIPotle source code. It returns scalar minimum value for array. It assumes all values in array are numeric.

=cut
	my $arrayref = shift or die;
	my @unsorted = @{ $arrayref };
	my @sorted = sort { $b <=> $a } @unsorted;
	my $max_index = $#sorted;
	return $sorted[$max_index];
}

sub ErrorFunction {

=head2 ErrorFunction(x)

This is a subroutine name from the original ChIPotle source code. It returns the integral of Gauss' standard error function for 'x'

=cut
	my $x = shift;
	my $result;
	# Returns the integral of Gauss' standard error function for x
	use constant kMaxLoop => 200;
	use constant kPi => 3.14159265358979;
	use constant kTiny => 1E-100;
	if ($x <= 2) {$result = _erf_low($x);}
	else{$result = _erf_high($x);}
	return $result;
	
	sub _erf_low {
		my $x = shift;
		my $t = 2 * $x * $x;
		my ($p, $s) = (1,1);
		for (my $i = 3; $i <= kMaxLoop; $i = $i + 2) {
			$p = $p * $t / $i;
			$s = $s + $p;
			if ($p < kTiny) { next }
		}
		my $y = 1 - 2 * $s * $x * exp((-1*$x) * $x) / sqrt(kPi);
		return $y;
	}
	
	sub _erf_high {
		my $x = shift;
		my ($a0, $b0, $a1, $b1, $f1) = (0, 0, 0, 1, 0);
		my ($a2, $b2, $f2);
		for (my $i = 1; $i <= kMaxLoop; $i = $i + 1) {
			my $g = 2 - ($i % 2);
			$a2 = $g * $x * $a1 + $i * $a0;
			$b2 = $g * $x * $b1 + $i * $b0;
			$f2 = $a2/ $b2;
			my $d = abs($f2 - $f1);
			if ($d < kTiny) { next };
			$a0 = $a1;
			$b0 = $b1;
			$a1 = $a2;
			$b1 = $b2;
			$f1 = $f2;
		}
		my $y = 2 * exp((-1*$x) * $x) / (2 * $x + $f2) / sqrt(kPi);
		return $y;
	}
}
