#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use Fasta;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use Log::Log4perl qw/:easy/;

Log::Log4perl->easy_init( { 
    level    => $DEBUG,
    #file     => ">run.log",
    layout   => '%d{HH:mm:ss} %p> (%L) %M - %m%n',
} );
my $logger = get_logger();

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
unless ($opt_reference_a && $opt_reference_b);

$logger->info("loading $opt_reference_a");
my $refa = slurp_fasta($opt_reference_a);

$logger->info("loading $opt_reference_b");
my $refb = slurp_fasta($opt_reference_b);

my @seqs = sort keys %$refa;
say "#Seq\tCoord\t$opt_reference_a\t$opt_reference_b";

my $counter = 0;

for my $seq (@seqs) {
    $logger->info("$counter") if ++$counter % 10000 == 0;
    $logger->info($seq);
    my $seqa = $refa->{$seq};
    my $seqb = $refb->{$seq};

    $logger->logdie("$seq lengths do not match") if length $seqa != length $seqb;
    my $length = length $seqa;

    $logger->info($length);

    for my $i (0..$length-1){
        my $aa = substr($seqa, $i, 1); 
        my $bb = substr($seqb, $i, 1);

        next if $aa eq $bb;

        printf "%s\t%d\t%s\t%s\n", $seq, $i + 1, $aa, $bb;
    }
}


=head1 NAME

genome_diff.pl - create diff between two genomes which 

=head1 SYNOPSIS

Usage examples:

 genome_diff.pl -a TAIR.fas -b Landsberg > diff

=head1 OPTIONS

=over

=item  -a <fasta> | --reference-a <fasta>

=for Euclid
    fasta.type: readable

=item  -b <fasta> | --reference-b <fasta>

=for Euclid
    fasta.type: readable

=back

=cut


