#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my $kmer = 2;
my @excluded_sequences;
my $alphabet = 'ACGT';
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'kmer|k=i'     => \$kmer,
    'exclude-seq|x=s{,}' => \@excluded_sequences,
    'alphabet|a=s' => \$alphabet,
    'output|o=s' => \$output,
    'verbose|v'  => sub { use diagnostics; },
    'quiet|q'    => sub { no warnings; },
    'help|h'     => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'   => sub { pod2usage ( -verbose => 2 ); }
);

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

my $reference = $ARGV[0];
my %reference = %{ index_fasta ($reference, join "|", @excluded_sequences) };
my %nucleotide_frequencies;

for my $k (1, $kmer) {
    
    print STDERR "#k:\t$k\n";

    my %frequencies  = ();
    my $total_length = 0;

    for my $chromosome (sort keys %reference) {
        next if $chromosome =~ m/-length/;

        $total_length += $reference{"$chromosome-length"};
        
        print STDERR "#chromosome:\t$chromosome\n";
        print STDERR "#length:\t", $reference{"$chromosome-length"}, "\n";

        my %frequency = %{ word_composition ($reference{$chromosome}, $k) };

        for my $word (sort keys %frequency) {

            next if $word =~ m/[^$alphabet]/i;

            $frequencies{$word} += $frequency{$word};

            print STDERR join ("\t",
                               $word,
                               $frequency{$word},
                               $frequency{$word} / $reference{"$chromosome-length"},
                           ), "\n";
        }
    }

    print STDERR "#all:\n";
    print STDERR "#length:$total_length\n";

    print join ("\t",
                '#word',
                '#count',
                '#size',
                '#observed',
                '#expected',
                '#obs/expect',
                '#independent',
            ), "\n"
            if $k > 1;

    for my $word (sort keys %frequencies) {

        next if $word =~ m/[^$alphabet]/i;

        if ($k == 1) {
            $nucleotide_frequencies{$word} = $frequencies{$word} / $total_length;
            print STDERR join ("\t",
                        $word,
                        $nucleotide_frequencies{$word},
                    ), "\n";
        }
        else {
            my $observed = $frequencies{$word} / $total_length;
            my $expected = 1;
            map {$expected *= $nucleotide_frequencies{$_}} (split //, $word);

            print join ("\t",
                        $word,
                        $frequencies{$word},
                        $total_length,
                        $observed,
                        $expected,
                        $observed / $expected,
                        map {$nucleotide_frequencies{$_} * $total_length} (split //, $word),
                    ), "\n";
        }
    }
}



sub word_composition {
    my ($sequence, $k, $alphabet) = @_;

    my %frequency = ();

    for (0 .. length ($sequence) - $k) {
     
        my $word = substr $sequence, $_, $k;
        next if $word =~ m/[^ACGT]/i;
        $frequency{$word}++;

    }

    return \%frequency;
}


sub index_fasta {
    my $reference_file     = shift;
    my $excluded_sequences = shift;

    my %reference = ();

    return \%reference unless $reference_file;

    # reads in the reference genome file into @fastaseq
    open my $REF, '<', "$reference_file"
    or croak "Can't open $reference for reading: $!";
    my @fastaseq = <$REF>;
    close $REF;

    # find and store indices for each chromosome change and corresponding descriptions
    my ( @idx, @dsc ) = ();
    for my $i ( 0 .. @fastaseq - 1 ) {
        if ( $fastaseq[$i] =~ m/^>/ ) {
            $fastaseq[$i] =~ s/>//g;
            $fastaseq[$i] = ( split /\s/, "$fastaseq[$i]" )[0];
            next if grep {m/$excluded_sequences/i} $fastaseq[$i];
            push @idx, $i;
            push @dsc, $fastaseq[$i];
        }
    }

    for my $j ( 0 .. @idx - 1 ) {
        my $line;
        if ( $j == scalar @idx - 1 ) {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. @fastaseq - 1]);
        }
        else {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. $idx[$j + 1] - 1]);
        }
        $line =~ s/[\n\r]//g;
        
        my $length = $line =~ tr/ACGT//;

        $reference{$dsc[$j]} = $line;
        $reference{"$dsc[$j]-length"} = $length;
    }
    return \%reference;
}

__END__


=head1 NAME

 name.pl - Short description

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

 name.pl [OPTION]... [FILE]...

 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: $:
 $Author: $:
 $Date: $:
 $HeadURL: $:
 $Id: $:

=head1 AUTHOR

 Pedro Silva <psilva@nature.berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
