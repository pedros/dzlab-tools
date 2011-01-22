#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

# Grabs and parses command line options
my $reference;
my @filter;
my $use_max_win = 0;
my $distance    = 0;
my $tss         = 0;
my $id          = 'ID';
my $use_strand;
my $gff;
my $default_locus; 
my $output;

my $result = GetOptions (
    'reference|r=s'     => \$reference,
    'filter|f=i{2}'     => \@filter,
    'use-max-win|w'     => \$use_max_win,
    'distance|d=i'      => \$distance,
    'tss|t'             => \$tss,
    'id|i=s'            => \$id,
    'default-locus|l=s' => \$default_locus,
    'use-strand|s'      => \$use_strand,
    'gff|g'             => \$gff,
    'output|o=s'        => \$output,
    'verbose|v'         => sub { use diagnostics; },
    'quiet|q'           => sub { no warnings; },
    'help|h'            => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'          => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result and $reference;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't read $output: $!";
    select $USER_OUT;
}

my %reference = %{ index_fasta ($reference) };

while (<>) {
    chomp; # delete line feeds
    next if ($_ =~ m/^#.*$|^\s*$|\.;\.$/);

    my %site = %{ gff_read ($_) };

    if ($use_max_win) {
        ($site{start}, $site{end})
        = $site{attribute} =~ m/
                                   maxstart = (\d+)
                                   .*
                                   maxend   = (\d+)
                               /xms;
    }

    my $length = $site{end} - $site{start};

    if ($distance) {
        my $center = $tss ? ( $site{strand} eq q{-} ? $site{end} : $site{start} )
                          : int ($length / 2) + $site{start};
        $site{start} = $center - int ($distance / 2);
        $site{start} = 0 if $site{start} < 0;
        $site{end}   = $center + int ($distance / 2);
        $length      = $distance;
    }

    if (@filter) {
        next if $site{attribute} eq q{.}
        or $length <= $filter[0]
        or $length >= $filter[1];
    }

    my ($attribute)
    = $site{attribute} =~ m/(?:\*|$id[\s=]?)([\w]+)/;

    $attribute ||= $default_locus // 'unknown_locus';

    {
        my ($seqname, $start) = @site{qw/seqname start/};
        my $sequence = substr $reference{$seqname}, ($start - 1), $length;

        if ($use_strand and q{-} eq $site{strand}) {
            $sequence = reverse $sequence;
            $sequence =~ tr/ACGTacgt/TGCAtgca/; 
            $sequence =~ s/\w$//;
        }
        elsif ($use_strand) {
            $sequence =~ s/^\w//;
        }

        if ($gff) {
            print join( "\t",
                        $site{seqname}, 'dz',       "$site{feature}_seq",
                        $site{start},   $site{end}, $site{score},
                        $site{strand},  q{.},       "seq=$sequence",
                    ), "\n";
        }
        else {
            print '>',$attribute
            ? "$attribute|"
            : '',
            "$site{seqname} $site{start}:$site{end}:$site{strand}\n";
            print "$sequence\n"
        }
    }
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
	'seqname'   => lc $seqname,
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

 # extract 100bp of sequence around transcription start site in each loci in annotation file
 ./extract_sequences.pl --reference path/to/genome.fasta annotation-file.gff --distance 100 --tss

=head1 DESCRIPTION

 Prints to STDOUT a fasta file with one sequence per line in the input gff file.
 The fasta header contains the locus id, and any other attributes, plus the sequence id and coordinates.

=head1 OPTIONS

 -f, --filter           only include sequences of length between -f n m
 -w, --use-max-win      try to find attribute fields 'maxstart=n' and 'maxend=m'
 -d, --distance         only go distance -d n/2 from center of locus 
 -i, --id               locus gff ID tag
 -t, --tss              extract $distance around transcription start site
 -r, --reference        reference fasta file from which to extract sequences.
 -l, --default-locus    default fasta header to use when extracting locus ID is not possible
 -s, --use-strands      reverse complement sequences as needed
 -g, --gff              format output in GFF format
 -o, --output           filename to write results to (defaults to STDOUT)
 -v, --verbose          output perl's diagnostic and warning messages
 -q, --quiet            supress perl's diagnostic and warning messages
 -h, --help             print this information
 -m, --manual           print the plain old documentation page 

=head1 VERSION

 $Rev: 395 $:
 $Author: psilva $:
 $Date: 2010-08-13 12:33:24 -0700 (Fri, 13 Aug 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/extract_sequences.pl $:
 $Id: extract_sequences.pl 395 2010-08-13 19:33:24Z psilva $:

=head1 REVISION

 0.0.2

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
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
