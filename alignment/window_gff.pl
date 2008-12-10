#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
disable diagnostics;

use Data::Dumper;
use Getopt::Long;
use File::Temp qw/ :mktemp  /;

# Globals, passed as command line options
my $gfffile = "";
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
    "gff|f=s" => \$gfffile,
    "width|w=i" => \$width,
    "step|s=i" => \$step,
    "output|o:s" => \$output,
    "verbose|v" => sub {enable diagnostics;use warnings;},
    "quiet|q" => sub {disable diagnostics;no warnings;},
    "usage|help|h" => \&usage
    );

# redirects STDOUT to file if specified by user
if(!($output eq '-')) {
    open(STDOUT, ">", "$output") or die("Can't redirect STDOUT to file: $output");
}

# opens gff file
open(my $GFF, "<", $gfffile) or die("Can't read file: $gfffile");

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

print STDERR "Sorting $gfffile by chromosomes, features, starting and ending coordinates";
my @data = <$GFF>;
@data = &gff_sort(@data);

print STDERR @data;exit;

print STDERR "Windowing data at $width bp width with $step bp intervals";
for(my $i = 0; $i < @data;) {

    my ($c_count, $t_count, $score) = (0, 0, 0);

    next if($data[$i] =~ m/\s*#/);

    my $start = ${&readGFF($data[$i])}{'start'};
    print STDERR "START: $start";
    my $end = $start + $width;
    my ($chromosome, $context);

    foreach my $line (@data[$i..($i + $width)]) {		
	chomp($line);
	next if($line =~ m/\s*#/);

	my %record = %{&readGFF($line)};

	last if( ($record{'seqname'} ne $chromosome) or
		 ($record{'feature'} ne $context) );

	$chromosome = $record{'seqname'};
	$context = $record{'feature'};

	my ($c_tmp, $t_tmp) = split(/;/, $record{'attribute'});
	$c_tmp =~ m/\d+/;
	$t_tmp =~ m/\d+/;
	$c_count += $c_tmp;
	$t_count += $t_tmp;
    }

    print join("\t",
	       $chromosome,
	       "avg",
	       ".",
	       $start,
	       $end,
	       $c_count / ($c_count + $t_count),
	       ".",
	       ".",
	       "c=$c_count;t=$t_count",
	), "\n";
    $i += $step;
}
print STDERR "Done";

close($GFF);
close(STDOUT);


# gff_sort sorts gff lines by sequence name, feature, start, and end coordinates
#-----------------sortArray-------------------------------
# sorts input array at fourth field delimited by tabs (starting coordinates)
# returns array sorted numerically on the 4th field
sub gff_sort {
    return sort {
	(split '\t', $a)[0] cmp (split '\t', $b)[0] or
	    (split '\t', $a)[2] cmp (split '\t', $b)[2] or
	    (split '\t', $a)[3] <=> (split '\t', $b)[3]
    } @_;
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
