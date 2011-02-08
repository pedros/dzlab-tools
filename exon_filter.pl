#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/all/;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use Pod::Usage;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/DZLab-Tools/lib";
use DZLab::Tools::GFF qw/gff_to_string gff_slurp/;

my $locus_tag = q/ID/;
my $input;
my $output_base;
my $help;
my $all;
my $result = GetOptions (
    "locus-tag|l=s"   => \$locus_tag,
    "output-base|o=s" => \$output_base,
    "input|i=s"       => \$input,
    "all"             => \$all,
    "help"            => \$help,
);
pod2usage(-verbose => 1) if (!$result || $help || !$input || !$output_base);  

# return hashref { locus_id => [gff_hashes] }

sub parse_exons{
    my ($opt) = @_;
    if (ref $opt->{gff_aref} ne 'ARRAY' || !exists $opt->{tag}){
        croak "parse_exon argument error";
    }
    my %exon;
    my $counter=0;
    my $no_locus_count=0;
    # populate %exon, { locus_id => [gff_hashrefs] }
    foreach my $gff (@{$opt->{gff_aref}}) {

        $counter++;
        if (! exists $gff->{$opt->{tag}}){
            $no_locus_count++;
            next;
        }

        # if ID=AT12345.123, use AT12345 part. else use entire thing
        # add an exonid entry into $gff
        my $locus;
        if ($gff->{$opt->{tag}} =~ /([^\.]+)\.([^\.]+)/){
            $locus = $1;
            $gff->{exonid} = $2;
        } else {
            $locus = $gff->{$opt->{tag}};
            carp "non-exon locus $locus already exists?"
            if exists $exon{$locus};
            $gff->{exonid} = '1';
        }
        push @{$exon{$locus}}, $gff;
    }
    carp "WARNING: $no_locus_count out of $counter gff records had no valid locus tag"
    if $no_locus_count;

    # sort, taking strand into account
    foreach my $locus (keys %exon) {
        my @strands = map { $_->{strand} // '+' } @{$exon{$locus}};
        my $dir = $strands[0];
        if (! all { $_ eq $dir } @strands){
            die "$locus: mixed strands for an exon?";
        }
        # sort backwards if - strand
        if ($dir eq '+'){
            @{$exon{$locus}} = sort { $a->{start} <=> $b->{start} } @{$exon{$locus}};
        } else {
            @{$exon{$locus}} = sort { $b->{start} <=> $a->{start} } @{$exon{$locus}};
        }
    }

    # if unique is set, only keep one version of the exon, chosen randomly.
    if ($opt->{unique}){
        foreach my $locus (keys %exon) {
            my @exonids = map { $_->{exonid} } @{$exon{$locus}};
            my $chosen_one = $exonids[0];
            @{$exon{$locus}} = grep { $_->{exonid} eq $chosen_one} @{$exon{$locus}};
        }
    }

    return \%exon;
}

my $gff_aref = gff_slurp({file => $input});

my $exon = parse_exons({gff_aref => $gff_aref,tag => $locus_tag,unique => ! $all});

open my $first  , '>' , "$output_base.first.gff";
open my $middle , '>' , "$output_base.middle.gff";
open my $last   , '>' , "$output_base.last.gff";

foreach my $locus (sort keys %$exon) {
    my $length = scalar @{$exon->{$locus}};
    say $first gff_to_string $exon->{$locus}[0];

    if ($length >= 3){
        foreach my $i (1 .. $length - 2) {
            say $middle gff_to_string $exon->{$locus}[$i];
        }
    }
    if ($length >= 2){
        say $last gff_to_string $exon->{$locus}[-1];
    }
}

close $first;
close $middle;
close $last;

=head1 NAME

exon_filter.pl - extract first, last, and middle exons for each locus into separate files

=head1 SYNOPSIS

Split TAIR8_exons.gff exons into TAIR8_exons.first.gff, TAIR8_exons.middle.gff and TAIR8_exons.last.gff

 exon_filter.pl -i TAIR8_exons.gff -o TAIR8_exons --all -l Parent

=head1 OPTIONS

 --help         print this information

 --all          output exons from all models for each loci.  default is 
  -a            to output only a single, arbitrary chosen model for each loci

 --locus-tag    Locus tag (default: 'ID').  
  -l  

 --input        name of input exon gff file (required)
  -i

 --output-base  base name of output files.  files will be named 
  -o            BASENAME.first.gff, BASENAME.middle.gff, and BASENAME.last.gff


=cut

