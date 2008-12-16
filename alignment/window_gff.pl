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
my @argv = @ARGV;

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
my @data;
while (<$GFF>) {
    chomp;
    next if ($_ =~ m/^#.*$|^\s*$/);
    push @data, $_;
}
close($GFF) if ($GFF ne "STDIN");

# prints out header fields that contain gff v3 header, generating program, time, and field names
&gff_print_header ($0, @argv);

# gets indices for chromosome changes (tracks changes on 0th field - seqname)
@data =  &gff_sort (@data);
my @locsindex = &gff_find_array ( 0, @data );

#print STDERR join("\n", @locsindex);exit;

for (my $i = 0; $i < @locsindex; $i++) {

    my @subdata;
    if ($i+1<@locsindex) {
	@subdata = @data[$locsindex[$i]..$locsindex[$i+1]];
    }
    else {
	@subdata = @data[$locsindex[$i]..scalar(@data) - 1];
    }

    # gets indices for context changes (tracks changes on 2th field - feature)
    my @featureindex = &gff_find_array ( 2, @subdata );

    for (my $l = 0; $l < @featureindex; $l++) {
	if ($l+1<@featureindex) {
	    &gff_sliding_window ( $width, $step, @subdata[$featureindex[$l]..$featureindex[$l+1]]);
	}
	else {
	    &gff_sliding_window ( $width, $step, @subdata[$featureindex[$l]..scalar(@subdata) - 1]);
	}
    }
}

close(STDOUT);


sub gff_sliding_window {
    my ($width, $step, @data) = @_;

    print STDERR
	"windowing data on sequence ${&gff_read ($data[scalar(@data)-1])}{'seqname'} and context ${&gff_read ($data[scalar(@data)-1])}{'feature'} with $width bp width and $step bp step...\n";

    my $lastcoord = ${&gff_read ($data[scalar(@data)-1])}{'end'};
    my $seqname = ${&gff_read ($data[scalar(@data)-1])}{'seqname'};
    my $context = ${&gff_read ($data[scalar(@data)-1])}{'feature'};
    my $lastrecord = 0;

    for (my $i = 1; $i < $lastcoord; $i++) {
	
	my ($c_count, $t_count, $score) = (0, 0, -0.1);

	my @range = &gff_filter_by_coord ($i, $i + $width, $lastrecord, \@data);

	$lastrecord = shift (@range);
	
	foreach my $k (@range) {
	    my %record = %{&gff_read ($k)};
	    my ($c_tmp, $t_tmp) = split(/;/, $record{'attribute'});
	    $c_tmp =~ m/(\d+)/;
	    $c_count += $1;
	    $t_tmp =~ m/(\d+)/;
	    $t_count += $1;
	    $seqname = $record{'seqname'};
	    $context = $record{'feature'};
	}

	if ($c_count + $t_count != 0) {
	    $score = $c_count / ($c_count + $t_count) unless $c_count == 0;
	}

	my $attribute = "c=$c_count;t=$t_count";

	if (scalar(@range) == 0) {
	    $attribute = ".";
	    $score = 0;
	}

	print join("\t",
		   $seqname,
		   "avg",
		   $context,
		   $i,
		   $i + $width - 1,
		   sprintf("%.3f", $score),
		   ".",
		   ".",
		   $attribute,
	    ), "\n";
	$i += $step - 1;
    }
}

# gff_sort sorts gff lines by sequence feature and start coordinate
sub gff_sort {
    print STDERR "sorting data...\n";
    return sort {
	(split '\t', $a)[0] cmp (split '\t', $b)[0] or
	(split '\t', $a)[2] cmp (split '\t', $b)[2] or
	(split '\t', $a)[3] <=> (split '\t', $b)[3]
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
    my ($lower, $upper, $last, $dataref) = @_;
    my (@data, @filtered) = @{$dataref};

    for (my $i = $last; $i < @data; $i++) {
	my %record = %{&gff_read ($data[$i])};
	if ($record{'start'} >= $lower && $record{'start'} <= $upper) {
	    push @filtered, $data[$i];
	    $last = $i;
	}
	last if ( $record{'start'} > $upper);
    }
    unshift (@filtered, $last);
    return @filtered;
}



#-----------Find indices for each ID change---------------#
# finds each index that signifies a change in type of record
# returns index array
sub gff_find_array {
    my ($field, @array) = @_;

    print STDERR "finding indices...\n";

    # @index contains a list of locations where array should be split
    # $previousid and $currentid are scalars containing the previous and current IDs (chr1, 2, etc)
    # $chrcount is just a counter for the number of different types of IDs
    my (@index, $previous, $current);
    my $chrcount = 0;
    
    # goes through full gene file
    for (my $i = 0; $i <@array ; $i++) {

	# gets current starting coordinate
	$current = (split '\t', $array[$i])[$field]; #chr1, chr2, etc

	# if we're at beginning of file if doesn't make sense to look for changes already
	# gets previous starting coordinate
	if ($i != 0) {$previous = (split '\t', $array[$i-1])[$field];}	
	else {$previous = $current;}
	
	# keeps track of number of different types of records
	# also stores each record type change in @index
	# ignores pound (#) characters if they're the first printing character in the record
	if ( ($current ne $previous) || ($i == 0) )  {
	    $index[$chrcount] = $i;
	    $chrcount++;
	}
    }
    return @index;
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
