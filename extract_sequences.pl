#!/usr/bin/perl

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
my $output;

my $result = GetOptions (
    'reference|r=s' => \$reference,
    'filter|f=i{2}' => \@filter,
    'use-max-win|w' => \$use_max_win,
    'distance|d=i'  => \$distance,
    'output|o=s'    => \$output,
    'verbose|v'     => sub { use diagnostics; },
    'quiet|q'       => sub { no warnings; },
    'help|h'        => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'      => sub { pod2usage ( -verbose => 2 ); }
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
        my $center = int ($length / 2) + $site{start};
        $site{start} = $center - ($distance / 2);
        $site{end}   = $center + ($distance / 2);
        $length      = $distance;
    }

    if (@filter) {
        next if $site{attribute} eq q{.}
        or $length <= $filter[0]
        or $length >= $filter[1];
    }

    my ($attribute)
    = $site{attribute} =~ m/\*([\w]+):/;

    $attribute ||= 'unknown_locus';

    print ">$attribute|$site{seqname}:$site{start}:$site{end}\n";
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

=head1 DESCRIPTION

 Prints to STDOUT a fasta file with one sequence per line in the input gff file.
 The fasta header contains the locus id, and any other attributes, plus the sequence id and coordinates.

=head1 OPTIONS

 -f, --filter      only include sequences of length between -f n m
 -w, --use-max-win try to find attribute fields 'maxstart=n' and 'maxend=m'
 -d  --distance     only go distance -d n / 2 from center of locus 
 -r, --reference   reference fasta file from which to extract sequences.
 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page 

=head1 VERSION

 $Rev$:
 $Author$:
 $Date$:
 $HeadURL$:
 $Id$:

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
