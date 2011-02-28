#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Log::Log4perl qw/:easy/;
use English;
Log::Log4perl->easy_init({ 
        level    => $DEBUG,
        layout   => "%d{HH:mm:ss} %p> (%L) %M - %m%n", 
        #file     => ">>log4perl.log",
    });

use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ! ($opt_input);

my $header_seq = substr $opt_adapter, 0, $opt_header_length;
my $header = q/$header_seq/;

DEBUG("$header_seq");

my ($counter, $too_short, $no_header) = (0,0,0);

open my $fh, '<', $opt_input;
my $output_file = $opt_out || $opt_input . '.trimmed';
open my $outfh, '>', $output_file;


while (defined (my $line = readfastq($fh))){
    INFO("$counter") if ++$counter % 100000 == 0;
    if ($line->[1] =~ m/$header_seq/){
        my $p = $PREMATCH;
        my $len = length $p;

        if ($len < $opt_minimum_insert_size){
            ++$too_short;
        } else{
            printf $outfh "%s\n%s\n%s\n%s\n", $line->[0], $p, $line->[2],  (substr $line->[3], 0, $len);
        }
    } else{
        ++$no_header;
    }
}
close $fh;
close $outfh;

INFO("total number of read: $counter");
INFO("inserts that were too short: $too_short");
INFO("inserts that did not have the adapter header: $no_header");

sub readfastq{
    my $fh = shift;
    my @line = map { my $x = scalar <$fh>; $x =~ tr/\r\n//d if defined $x; $x } (1..4);
    if (grep { defined } @line){
        return \@line;
    }
    return;
}

=head1 NAME

fastq_trim.pl - Trim out adapter for fastq reads.

=head1 SYNOPSIS

trim off adapter by finding 10 bp header from default adapter sequence and cutting it and trailing bases off:

 fastq_trim.pl -i input.fastq -h 10 

same as above but use a custom adapter sequence:

 fastq_trim.pl -i input.fastq -h 10 -a ACCCGGTGAGACATGAC

Even when you pass a custom adapter, only '-h' bases will be used.

=head1 OPTIONS

=over

=item  --adapter <seq> | -a <seq>

Adapter sequence.  By default it is:

 AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG

=for Euclid
    seq.default: "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG"

=item --header-length <len> | -h <len>

Number of bases from the adapter sequence to use to search.  Default is 7.

=for Euclid
    len.type: int > 0 
    len.default: 7

=item  --minimum-insert-size <size> | -m <size>

Don't report files under this long. Default to 39.

=for Euclid
    size.type:    int > 0
    size.default: 39

=item --input <file> | -i <file>

Input file

=item --out <file> | -o <file>

Output file.  Default to <input_file>.trimmed

=back

=cut

