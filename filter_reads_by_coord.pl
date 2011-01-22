#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

# print usage if no cmd parameters
pod2usage ( -verbose => 1 )
unless @ARGV > 0;

my $gff_filter; # gff file with coordinates to filter by
my $gff_reads;  # gff file with reads from correlatePairedEnds.pl
my $reference;  # fasta file for calculating chromosome lengths for reversing coords
my $output;     # output file name

# Grabs and parses command line options
my $result = GetOptions (
    'gff-filter|f=s' => \$gff_filter,
    'gff-reads|g=s'  => \$gff_reads,
    'reference|r=s'  => \$reference,
    'output|o=s'     => \$output,
    'verbose|v' => sub { use diagnostics; },
    'quiet|q'   => sub { no warnings; },
    'help|h'    => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'  => sub { pod2usage ( -verbose => 2 ); }
);

# print usage if required parameters are missing
pod2usage ( -verbose => 1 )
unless $gff_filter and $gff_reads and $reference;

# hash of ints 'reference' holds all chromosome lengths
# hash of array refs 'filter' holds all coordinates to filter by
my %reference = %{ index_fasta ($reference) };
my %filter = ();

# index the gff filtering file
# optimized for specific gff file of with 10 fields (9th for c+t count)
open my $GFF_FILTER, '<', $gff_filter or croak "Can't read $gff_filter.";
FILTER:
while (<$GFF_FILTER>) {
    next FILTER if m/^\s*#/;

    my ($chr, $start, $end, $strand, $ct) =
    (split /\t/, $_)[0, 3, 4, 6, 8];

    next if $ct < 20; # disregard low coverage ranges

    # this is supposed to properly convert reverse strand coords
    # however, it is not working properly at this time, so I'm skipping reverse strand reads
    # $start = $reference{$chr} - $start + 1 if $strand eq q{-};
    # $end   = $reference{$chr} - $end   + 1 if $strand eq q{-};

    # if we're dealing with a range and there's strand information save the strand
    if ($start != $end and $strand =~ m/\+|-/) {
        $filter{$chr}{$start} = [$end, $strand];
    }
    # for ranges of coordinates, index all inner coordinates as keys (for speed)
    else {
        $filter{$chr}{$_} = [$end] for ($start .. $end);
    }

}

# start processing reads
open my $GFF_READS, '<', $gff_reads or croak "Can't read $gff_reads.";
my ($CG_READS, $TG_READS);
# if we're not output to STDOUT, open two streams, one for methylated reads, one for unmethylated reads
unless ($output eq q{-}) {
    open $CG_READS, '>>', "$output-cg-filter.gff" or croak "Can't write $output-cg-filter.gff";
    open $TG_READS, '>>', "$output-tg-filter.gff" or croak "Can't write $output-tg-filter.gff";
}
else {
    open STDOUT, '>', '-' or croak "Can't redirect to STDOUT.";
}
READ:
while (my $read = <$GFF_READS>) {
    # skip gff comments
    next READ if $read =~ m/^\s*#/;

    # match is number of alignment targets for this read
    my ($chr, $start, $end, $match, $strand) =
    (split /\t/, $read)[0, 3, 4, 5, 6];

    # only accept unique matches
    next READ unless $match == 1 and $strand eq '+'; # disregarding reverse strand reads due to problem processing those --> needs fixing

    # loop through each coordinate position in read
    # find whether current read overlaps one of the indexed coordinates
    POSITION:
    for my $pos ($start .. $end) {

        # if current position was indexed from the filtering file
        if (exists $filter{$chr}{$pos}) {

            # calculate overlap
            # if read terminates before indexed range end overlap is range end minus current pos
            # if read terminates after  indexed range end overlap is read  end minus current pos
            # skip reads that don't overlap at least by 70%
            my $overlap = 0;
            if ($end <= $filter{$chr}{$pos}->[0]) {$overlap = abs ($end - $pos)}
            else {$overlap = abs ($filter{$chr}{$pos}->[0] - $pos)}
            next READ if $overlap < abs ($end - $start) * 0.7;

            # if there is strand information available, skip read if not in same strand
            next POSITION
            if defined $filter{$chr}{$pos}->[1] and $filter{$chr}{$pos}->[1] ne $strand;

            # extract scaffold and read sequences. dash is /1 or /2 (direction of sequencing)
            my ($scaffold) = $read =~ m/target=([A-Z]+)/;
            my ($dash, $methread) = $read =~ m/\/([12]):([A-Z]+)/;

            # output selection method, if not printing to STDOUT (in which case we don't split)
            unless ($output eq q{-}) {

                my @methylated   = split //, $methread;
                my @unmethylated = split //, $scaffold;

                READ_COORD:
                for (my $i = $start; $i < $end; ++$i) {

                    # sets $j to 0 for forward strand coordinates
                    my $j = $i - $start;

                    # Since the input coordinates are absolute (ie. from the forward strand)
                    # if a particular read maps to the reverse strand, the coordinates for any
                    # 'c's hit need to be reversed
                    my $coord = $i;
                    $coord = $reference{$chr} - $i + 1 if $strand eq '-';

                    # checks that we're looking at a left/1 sequence AND current character is a 'C'
                    # this regex is adapted to the Solexa sequences we have, so it might change in the future
                    # we're looking 2 characters ahead in the scaffold because it has 2 extra chars in the beginning
                    if ( $unmethylated[$j + 2] =~ m/[Cc]/
                         and $dash == 1) {

                        # checks the context by looking ahead +3 or +3..+4 bps
                        # because the scaffold read is displaced by 2 bps
                        if ( $unmethylated[$j+3] =~ m/[Gg]/ ) {

                            # checks what happened in the methylated sequence
                            # and updates the appropriate 'c' or 't' count
                            if ( $methylated[$j] =~ m/[Cc]/ ) {
                                select $CG_READS;
                                last READ_COORD;
                            }
                            elsif ($methylated[$j] =~ m/[Tt]/ ) {
                                select $TG_READS;
                                next READ_COORD;
                            }
                        }
                    }

                    # checks that we're looking at a right/2 sequence AND current character is a 'C'
                    # this regex is adapted to the Solexa sequences we have, so it might change in the future
                    # we're looking 2 characters before in the scaffold because it has 2 extra chars in the beginning
                    # AND we're looking in the reverse strand, so the context goes right-to-left
                    elsif ( $unmethylated[$j + 2] =~ m/[Gg]/
                            and $dash == 2) {

                        # checks the context by looking behind +1 or +0..+1 bps
                        # because the scaffold read is displaced by 2 bps
                        if ( $unmethylated[$j + 1] =~ m/[Cc]/ ) {

                            # checks what happened in the methylated sequence
                            # and updates the appropriate 'c' or 't' count
                            if ( $methylated[$j] =~ m/[Gg]/ ) {
                                select $CG_READS;
                                last READ_COORD;
                            }
                            elsif ( $methylated[$j] =~ m/[Aa]/ ) {
                                select $TG_READS;
                                next READ_COORD;
                            }
                        }
                    }
                }
            }
            # print read if there is a CG site in scaffold
            print $read if $scaffold =~ m/CG/;

            last POSITION; # if we got here, this read belongs, and we can skip the rest of the bases
        }
    }
}
# close all files
close $GFF_READS;
close $CG_READS if $CG_READS;
close $TG_READS if $TG_READS;

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

    # tries to find each chromosome's centrometer center coordinate
    for my $j ( 0 .. @idx - 1 ) {
        my $line;
        if ( $j == scalar @idx - 1 ) {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. @fastaseq - 1]);
        }
        else {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. $idx[$j + 1] - 1]);
        }
        $line =~ s/[\n\r]//g;
        $reference{$dsc[$j]} = length $line;
    }
    return \%reference;
}



__END__


=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

=head1 REVISION

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
