#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

use FindBin;
use lib "$FindBin::Bin/lib";
use Fasta;

use Pod::Usage;
use Getopt::Long;

my $help;
my $reference;
my $result = GetOptions (
    "reference|r=s" => \$reference,
    "help"    => \$help,
);
pod2usage(-verbose => 1) if (!$result || $help);  

my $genome = slurp_fasta($reference);

my $counter = 0;
while (defined(my $line = <ARGV>)){
    $line =~ tr/\r\n//d;
    my ($seq, $coord, $attr) = (split /\t/, $line)[0,3,8];
    
    next unless $seq && $coord && $attr;
    my ($original, $converted) = split />/,$attr;
    die "$line fubar" unless $original && $converted;

    $seq = uc $seq;

    # check before
    {
        my $grab = substr $genome->{$seq}, $coord-1,1;
        if ($grab ne $original){
            say STDERR "BEFORE: $grab doesn't match $original at $line"; 
            ++$counter;
            next ;
        }
    }

    substr($genome->{$seq}, $coord-1,1) = $converted;

    # check after
    {
        my $grab = substr $genome->{$seq}, $coord-1,1;
        warn "AFTER: $grab doesn't match $original at $line" unless $grab eq $converted;
    }
}

for my $seq (sort keys %$genome) {
    say format_fasta($seq, $genome->{$seq});
}
