#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use File::Basename;

use FindBin;
use lib "$FindBin::Bin/DZLab-Tools/lib";
use DZLab::Tools::Fasta;

use Pod::Usage;
use Getopt::Long;

my $help;
my $verbose;
my $input;
my $output=q/-/;
my $window_size = 50;
my $result = GetOptions (
    "input|i=s"       => \$input,
    "output|o=s"      => \$output,
    "window-size|w=i" => \$window_size,
    "help"            => \$help,
);
pod2usage(-verbose => 1) if (!$result || $help || !$window_size || !$input);  

my $sequences = slurp_fasta($input, {-l => 1});

my $outfh;
if ($output eq q{-}){
    $outfh = *STDOUT;
}
else {
    open $outfh, '>', $output;
}

foreach my $seqname (sort keys %$sequences) {
    my $seq = $sequences->{$seqname};
    my $length = length $seq;

    my $pos = 0; # zero indexed

    while ($pos < $length){
        my $window = substr $seq, $pos, $window_size;
        my %count = (
            a => ($window =~ tr/a//),
            c => ($window =~ tr/c//),
            g => ($window =~ tr/g//),
            t => ($window =~ tr/t//),
        );
        my $score = $count{c};
        my $attr = join ";", map { "$_=$count{$_}" } qw/a c g t/;
        
        printf $outfh "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
        $seqname,
        basename($input),
        "window",
        # b/c people want 1-indexed 
        $pos+1, 
        # no -1 here b/c the start position also counts as one...
        $pos+$window_size > $length ? $length : $pos+$window_size, 
        $score,
        q{+},
        q{.},
        $attr;
        $pos += $window_size;
    }
}



if ($output ne q{-}){
    close $outfh;
}
=head1 NAME

window_genome.pl - count number of bases in a genome by fixed-size windows,
output a gff line for each window with counts in the atttributes and 
score set to the c count

=head1 SYNOPSIS

window_genome.pl -w 50 -i TAIR_reference.fas -o output.txt

=head1 OPTIONS

 --window-size   -w  Size of window.
 --input         -i  Input file
 --output        -o  Output file (default to screen)
 --help              print this information

=cut


