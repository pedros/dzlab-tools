#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Pod::Usage;
use Getopt::Long;

use FindBin;
use lib "$FindBin::Bin/DZLab-Tools/lib";
use DZLab::Tools::GFF qw/gff_make_iterator/;

my $help;
my $verbose;
my $column = 'ID';
my $output = q{-};
my $input;
my $unique;
my $sort;
my $result = GetOptions (
    "column|c=s" => \$column,
    "output|o=s" => \$output,
    "input|i=s"  => \$input,
    "unique|u"   => \$unique,
    "sort|s"     => \$sort,
    "verbose"    => \$verbose,
    "help"       => \$help,
);
pod2usage(-verbose => 99) if (!$result || !$input || $help);  


# map numbered columns to named

my %default_cols = (
    1 => 'seqname',
    2 => 'source',
    3 => 'feature',
    4 => 'start',
    5 => 'end',
    6 => 'score',
    7 => 'strand',
    8 => 'frame',
    9 => 'attribute',
);
$column = $column =~ /^\d$/ ? $default_cols{$column} : $column;

unless ($output eq '-'){
    open STDOUT, '>', $output;
}

my $iter = gff_make_iterator(file => $input);

my @accum;
while (defined(my $row = $iter->())){
    my $val = $row->{$column};
    if (defined $val){
        push @accum, $val;
    }
}

if ($unique){
    my %uniq = map {$_ => 0} @accum;
    @accum = keys %uniq;
}
if ($sort) {
    @accum = sort @accum;
}

foreach my $row (@accum) {
    say $row;
}

=head1 NAME

extract_column_gff.pl - extract a single column or attribute field from gff

=head1 SYNOPSIS

This script can extract either entire columns (which can be denoted by their number 1 through 9, or by name 
seqname, source, feature, start, end, score, strand, frame, or attribute) or specific fields in the attributes field. 
For example, for the GFF line below, you can specify a column as ID, Name, or Note: 

 Chr1 TAIR8 gene 6790 8737 . - . ID=AT1G01020;Name=AT1G01020;Note=ARV1

Examples:

Grab all ID's from input.gff, sort it, get rid of duplicates, and put it into output:

 extract_column_gff.pl -i input.gff -o output.txt -s -u 

Grab sequence names (column 1) from input.gff, get rid of duplicates, print to screen:

 extract_column_gff.pl -i input.gff -c 1 -u 

=head1 OPTIONS

    --verbose      print increasingly verbose error messages
    --help         print this information

    --column  -c   column/field name to extract. (Default: ID)
                   can be a number 1 through 9, attribute field key name 
                   (such as 'ID', 'c', 'Note', etc).

    --output  -o   Output file. Defaults to stdin (screen).
    --input   -i   GFF File to filter. required.
    --unique  -u   Get rid of duplicates.
    --sort    -s   Sort alphabetically 
=cut

