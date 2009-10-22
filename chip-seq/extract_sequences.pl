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

# Grabs and parses command line options
my $reference;
my $filter;
my $result = GetOptions (
    'reference|r=s' => \$reference,
    'filter|f'      => \$filter,
    'verbose|v'     => sub { use diagnostics; },
    'quiet|q'       => sub { no warnings; },
    'help|h'        => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'      => sub { pod2usage ( -verbose => 2 ); }
);

my %reference = %{ index_fasta ($reference) };

while (<>) {
    chomp; # delete line feeds
    next if ($_ =~ m/^#.*$|^\s*$|\.;\.$/);

    my %site = %{ gff_read ($_) };

    my $length = $site{end} - $site{start};

    # next if $filter
    # and $site{attribute} eq q{.}
    # or $length <= 50
    # or $length >= 300;

    my ($attribute)
    = $site{attribute} =~ m/ID=([^;]+);/;
    $attribute =~ s/:/|/;

    print ">$attribute $site{seqname}:$site{start}:$site{end}\n";
    print substr $reference{$site{seqname}}, ($site{start} - 1), ($length);
    print "\n";
}


# index_fasta reads a standard fasta file and returns a hash with chromosome names as keys as chromosome lengths as values
sub index_fasta {
    my $reference_file = shift;

    # holds name of chromosomes as keys and length of chromosomes in bp as values
    my %reference = ();

    return \%reference unless $reference_file;

    # reads in the reference genome file into @fastaseq
    open my $REF, '<', "$reference_file" or croak "Can't open file: $reference_file";
    my @fastaseq = <$REF>;
    close $REF;

    # find and store indices for each chromosome change and corresponding descriptions
    my ( @idx, @dsc ) = ();
    for my $i ( 0 .. @fastaseq - 1 ) {
        if ( $fastaseq[$i] =~ m/^>/ ) {
            $fastaseq[$i] =~ s/>//g;
            $fastaseq[$i] = ( split /\s/, "$fastaseq[$i]" )[0];
            $fastaseq[$i] =~ tr/A-Z/a-z/;
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
        $reference{$dsc[$j]} = $line;
    }
    return \%reference;
}


sub gff_read {
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, shift);

    $seqname =~ tr/A-Z/a-z/;

    my %rec = (
	'seqname'   => $seqname,
	'source'    => $source,
	'feature'   => $feature,
	'start'     => $start,
	'end'       => $end,
	'score'     => $score,
	'strand'    => $strand,
	'frame'     => $strand,
	'attribute' => $attribute
	);
    return \%rec;
}


__END__


=head1 NAME

 extract_sequences.pl -- Given a gff file with start and end coordinates per line, extracts corresponding sequences from genome fasta file.

=head1 SYNOPSIS

 ./extract_sequences.pl --reference path/to/genome.fasta annotation-file.gff

=head1 DESCRIPTION

 Prints to STDOUT a fasta file with one sequence per line in the input gff file.
 The fasta header contains the locus id, and any other attributes, plus the sequence id and coordinates.

=head1 OPTIONS

 --reference, -r    Reference fasta file from which to extract sequences.

=head1 REVISION

 0.0.1

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
