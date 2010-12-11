#!/usr/bin/env perl

use strict;
use warnings;
use Test::Simple tests => 1;
use DZLab::Tools::Fasta;
use Test::Deep;
#use Test::More qw(no_plan);

my $width = 80;

sub write_random_fasta{
    my $file = shift || die;
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


my $tmpfile1 = "/tmp/tmp1.fasta";
my $tmpfile2 = "/tmp/tmp2.fasta";
write_random_fasta $tmpfile1;

my $read1 = slurp_fasta($tmpfile1);

open my $fh, '>', $tmpfile2;

while (my ($header,$seq) = each %$read1) {
    print $fh format_fasta($header,$seq);
}
close $fh;

my $read2 = slurp_fasta($tmpfile2);

cmp_deeply($read1, $read2, 
    "Equality of original fasta read and written and reread");

unlink $tmpfile1;
unlink $tmpfile2;
