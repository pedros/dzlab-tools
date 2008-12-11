#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;disable diagnostics;
use Data::Dumper;
use Getopt::Long;
use File::Temp qw/ :mktemp  /;

# Globals, passed as command line options
my $gfffile = "-";
my $width = 0;
my $step = 0;
my $output = "-";
my $verbose = 0;
my $quiet = 0;
my $usage = 0;

# Initial check of command line parameters
&usage;

# Grabs and parses command line options
my $result = GetOptions (
    "gff|f:s" => \$gfffile,
    "width|w=i" => \$width,
    "step|s=i" => \$step,
    "output|o:s" => \$output,
    "verbose|v" => sub {enable diagnostics;use warnings;},
    "quiet|q" => sub {disable diagnostics;no warnings;},
    "usage|help|h" => \&usage
    );

# redirects STDOUT to file if specified by user
if (!($output eq '-')) {
    open(STDOUT, ">", "$output") or die("Can't redirect STDOUT to file: $output");
}

# opens gff file or STDIN
my $GFF;
if (!($gfffile eq '-')) {open($GFF, "<", $gfffile) or die("Can't read file: $gfffile");}
else {$GFF = "STDIN";}

# reads in data
my @data = <$GFF>;
close($GFF) if ($GFF ne "STDIN");

# prints out header fields that contain gff v3 header, generating program, time, and field names
&gff_print_header ($0);

# gets indices for chromosome changes (tracks changes on 0th field - seqname)
my @locsindex = &gff_find_array ( 0, &gff_sort (@data) );

for (my $i = 0; $i < @locsindex; $i++) {

    my @subdata;
    if ($i+1<@locsindex) {
	 @subdata = &array_split ($locsindex[$i], $locsindex[$i+1], @data);
    }
    else {
	 @subdata = &array_split ($locsindex[$i], scalar(@data), @data);
    }

    # gets indices for context changes (tracks changes on 2th field - feature)
    my @featureindex = &gff_find_array ( 2, @subdata );

    for (my $l = 0; $l < @featureindex; $l++) {
	if ($l+1<@featureindex) {
	    	&gff_sliding_window ( $width, $step, &array_split ($featureindex[$l], $featureindex[$l+1], @subdata) );
	}
	else {
	    	&gff_sliding_window ( $width, $step, &array_split ($featureindex[$l], scalar(@subdata), @subdata) );
	}
    }
}

close(STDOUT);


sub gff_sliding_window {
    my ($width, $step, @data) = @_;

    my $lastcoord = ${&gff_read ($data[scalar(@data)-1])}{'end'};

    for (my $i = 1; $i < $lastcoord; $i++) {
	
	my ($c_count, $t_count, $score) = (0, 0, 0);
	
	my @range = &gff_filter_by_coord ($i, $i + $width, @data);
	
	foreach my $k (@range) {
	    my %record = &gff_read ($k);
	    my ($c_tmp, $t_tmp) = split(/;/, $record{'attribute'});
	    $c_tmp =~ m/\d+/;
	    $t_tmp =~ m/\d+/;
	    $c_count += $c_tmp;
	    $t_count += $t_tmp;
	}

	if ($c_count + $t_count != 0) {
	    $score = ( $c_count / ($c_count + $t_count) ) / $width;
	
	}

	print join("\t",
		   ${&gff_read ($data[$i])}{'seqname'},
		   "avg",
		   ${&gff_read ($data[$i])}{'feature'},
		   $i,
		   $i + $width,
		   $score,
		   ".",
		   ".",
		   "c=$c_count;t=$t_count",
	    ), "\n";
	$i += $step - 1;
    }
}

# gff_sort sorts gff lines by sequence feature and start coordinate
sub gff_sort {
    return sort {
	(split '\t', $a)[2] cmp (split '\t', $b)[2] or (split '\t', $a)[3] <=> (split '\t', $b)[3]
    } @_;
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
	'frame'=>$strand,
	'attribute'=>$attribute
	);
    return \%rec;
}


sub gff_filter_by_coord {

    my ($lower, $upper, @data) = @_;

    my @filtered;
    foreach my $i (@data) {
	my %record = %{&gff_read ($i)};
	push @filtered, $i if ( $record{'start'} >= $lower);
	last if ( $record{'start'} > $upper);
    }
    return @filtered;
}



#-----------Find indices for each ID change---------------#
# finds each index that signifies a change in type of record
# returns index array
sub gff_find_array {
    my ($field, @array) = @_;

    # @index contains a list of locations where array should be split
    # $previousid and $currentid are scalars containing the previous and current IDs (chr1, 2, etc)
    # $chrcount is just a counter for the number of different types of IDs
    my (@index, $previousid, $currentid);
    my $chrcount=0;
    
    # goes through full gene file
    for (my $k=0;$k<@array;$k++) {

	# gets current starting coordinate
	$currentid=(split '\t', $array[$k])[$field]; #chr1, chr2, etc
	
	# if we're at beginning of file if doesn't make sense to look for changes already
	# gets previous starting coordinate
	if ($k != 0) {$previousid = (split '\t', $array[$k-1])[$field];}	
	else {$previousid = $currentid;}
	
	# keeps track of number of different types of records
	# also stores each record type change in @index
	# ignores pound (#) characters if they're the first printing character in the record
	if ( ($currentid ne $previousid && $currentid !~ m/^\s*#/) || $k==0 )  {
	    if($currentid !~ m/^\s*#/) {
		$index[$chrcount]=$k;
		$chrcount++;
	    }
	}
    }
    return @index;
}

#----------Split into multiple arrays----------------#
# takes an input index array
# returns split array
sub array_split {
    my ($start, $end)=($_[0], $_[1]);
    my @array=@_[2..@_];
    return @array[$start..$end-1];
}

# prints out usage information
sub usage {
    if ($usage || @ARGV<2) {
	print STDERR
	    "countMethylation.pl <PARAMETERS> [OPTIONS]
\t<--gff>\tGFF alignment input file
\t<--width>\tWidth size of sliding window in bp
\t<--step>\tStep interval of sliding window in bp
\t[--output]\tFilename to write results to (default is STDOUT)
\t[--verbose]\tOutput perl's diagnostic and warning messages
\t[--quiet]\tSupress perl's diagnostic and warning messages
\t[--usage]\tPrint this information
Takes a GFF-formatted input file with methylation information.
Each line corresponds to a single 'c' in the genome that was sequenced.
Assumes that:
    1) input GFF file is sorted by starting coordinate AND by context
    2) input GFF file contains a single sequence id (ie. single chromosome)
Runs a sliding window of width x and step interval y and averages the score
for each step, generating a GFF file with N/y lines (N=number of input 'c's).
";
	exit 1;
    }
}


# prints out a commented line with header fields
sub gff_print_header {
    print "##gff-version 3\n";
    print "# ", shift, ": ";
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime (time);
    printf "%4d-%02d-%02d %02d:%02d:%02d\n", $year+1900, $mon+1, $mday, $hour, $min, $sec;
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
