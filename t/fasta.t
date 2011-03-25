#!/usr/bin/env perl
use Data::Dumper;
use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib";
use Fasta;

use Test::More qw(no_plan);
use feature 'say';


my $width = 80;

sub write_random_fasta{
    my $file = shift || die "need filename";
    open my $fh, '>', $file or die "can't open file";
    for (1 .. 100){
        print $fh ">seq$_ " . random_seq(12) . "\n";
        for (1.. 4 + int rand(100)){
            print $fh random_seq() . "\n";
        }
        print $fh random_seq(int rand 70) . "\n";
    }
    close $file;
}

sub random_seq{
    my $length = shift // $width;

    my @bases = qw/A C T G/;

    return join q{}, map { $bases[int rand(4)] } (1 .. $length); 
}

use feature 'switch';

sub check_pattern{
    my ($seq1, $seq2, $pattern) = @_;
    given ($pattern){
        when (q/c2t/){ $seq1 =~ tr/cC/tT/; }
        when (q/g2a/){ $seq1 =~ tr/gG/aA/; }
        when (q/rc/){ $seq1 =~ tr/acgtACGT/tgcaTGCA/; }
        default {die "unknown pattern"}
    }
    return $seq1 eq $seq2;
}

my $tmpfile1 = "/tmp/tmp1.fasta";
my $tmpfile2 = "/tmp/tmp2.fasta";
unlink $tmpfile1;
unlink $tmpfile2;

write_random_fasta $tmpfile1;

my $read1 = slurp_fasta($tmpfile1);

open my $fh, '>', $tmpfile2 or die "can't open temporary outfile";

while (my ($header,$seq) = each %$read1) {
    print $fh format_fasta($header,$seq);
}
close $fh;

my $read2 = slurp_fasta($tmpfile2);

is_deeply($read1, $read2, 
    "Equality of original fasta read and written and reread");

unlink $tmpfile1;
unlink $tmpfile2;
