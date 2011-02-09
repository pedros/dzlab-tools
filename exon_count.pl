#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/all/;
use List::Util qw/max min/;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use Pod::Usage;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;

my $locus_tag = q/ID/;
my $input;
my $output = q{-};
my $help;
my $result = GetOptions (
    "locus-tag|l=s" => \$locus_tag,
    "input|i=s"     => \$input,
    "output|o=s"    => \$output,
    "help"          => \$help,
);
pod2usage(-verbose => 1) if (!$result || $help || !$input || !$output);  

unless ($output eq '-'){
    close STDOUT;
    open STDOUT, '>', $output or die "can't open $output for writing";
}

my $p = GFF::Parser->new(file => $input,locus => $locus_tag);
my $index = $p->slurp_index($locus_tag);

while (my ($locus,$seqs) = each %$index) {
    my @exons = grep {$_->{feature} eq 'exon'} @$seqs;
    my $min = (min sort { $a->{start} <=> $b->{start} } @exons)->{start};
    my $max = (max sort { $a->{end} <=> $b->{end} } @exons)->{end};
    say join "\t", 
    $exons[0]->{seqname} || '.',
    'dzlab',
    'exon',
    $min,
    $max,
    scalar @$seqs,
    $exons[0]->{strand} || '.',
    '.',
    join("=",$locus_tag,$locus);
}

=head1 NAME

exon_count.pl - count number of exons for each gene model

=head1 SYNOPSIS

 exon_count.pl -l Parent -i TAIR8_gmod-exons.gff -o

=head1 DESCRIPTION

Long description

=head1 OPTIONS

 --locus-tag  Locus tag (default: 'ID').  
  -l  

 --input      name of input exon gff file (required)
  -i

 --output     output file. defaults to stdout (screen)
  -o          

 --help       print this information

=cut

