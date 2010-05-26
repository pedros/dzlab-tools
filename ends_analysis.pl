#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(min max sum);

my $DATA_HANDLE    = 'ARGV';
my $gff_annotation = q{};
my $bin_width      = 100;
my $distance       = 5000;
my $stop_flag      = 2;
my $stop_distance  = 1500;
my $three_prime    = 0;
my $five_prime     = 0;
my $attribute_id   = 'ID';
my $output;

# Grabs and parses command line options
my $result = GetOptions(
    'gff-annotation|g=s' => \$gff_annotation,
    'bin-width|b=i'      => \$bin_width,
    'distance|d=i'       => \$distance,
    'stop-flag|s=i'      => \$stop_flag,
    'stop-distance|k=i'  => \$stop_distance,
    'three-prime|3'      => \$three_prime,
    'five-prime|5'       => \$five_prime,
    'extract-id|x=s'     => \$attribute_id,
    'output|o=s'         => \$output,
    'verbose|v'          => sub { use diagnostics; },
    'quiet|q'            => sub { no warnings; },
    'help|h'             => sub { pod2usage( -verbose => 1 ); },
    'manual|m'           => sub { pod2usage( -verbose => 2 ); }
);

# Add additional stop strategies here and implement the corresponding function
my $stop_flag_dispatch = {
    0 => \&stop_flag_0,    # implemented; passed tests
    1 => \&stop_flag_1,
    2 => \&stop_flag_2,    # implemented; check embedded loci
    3 => \&stop_flag_3,
    4 => \&stop_flag_4,
    5 => \&stop_flag_5,
    6 => \&stop_flag_6,    # implemented: see #2
};

# Check required command line parameters
pod2usage( -verbose => 1 )
    unless @ARGV
    and $result
    and $gff_annotation
    and ( $three_prime xor $five_prime )
    and ( $stop_flag != 6 or ( $stop_flag == 6 and $stop_distance ) )
    and exists $stop_flag_dispatch->{$stop_flag},
    and $bin_width > 0;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't read $output: $!";
    select $USER_OUT;
}

# index and offset the loci annotations data based on stop flag
my $annotation = offset_gff_annotation(
    index_gff_annotation( $gff_annotation, $attribute_id ),
    $stop_flag_dispatch->{$stop_flag},
    {   -three_prime   => $three_prime,
        -five_prime    => $five_prime,
        -distance      => $distance,
        -stop_flag     => $stop_flag,
        -stop_distance => $stop_distance
    }
);

my $num_bins           = 2 * int( $distance / $bin_width );
my %genes              = (); # hash of loci id hashes with scores/strand keys
my %sorted_annotations = (); # cache for sorted GFF annotations
my $gff_iterator       = make_gff_iterator( $ARGV[0], \&gff_read );

COORD:
while ( my $gff_line = $gff_iterator->() ) {

    # Step 1: Parsing GFF line
    # &gff_read returns [] for GFF comments, invalid lines, etc.
    next COORD unless ref $gff_line eq 'HASH';

    next COORD if $gff_line->{score} eq q{.};

    # Step 2: Memoization
    # has this been done before? use it : otherwise cache it
    unless ( exists $sorted_annotations{ $gff_line->{seqname} } ) {

        my @loci_coords   = keys %{ $annotation->{ $gff_line->{seqname} } };
        my @sorted_coords = sort { $a <=> $b } @loci_coords;

        # map to values [start, end, strand, locus id]
        $sorted_annotations{ $gff_line->{seqname} } = [
            map { $annotation->{ $gff_line->{seqname} }{$_} }
                @sorted_coords
        ];
    }

    # Step 3: Searching loci where this score fits into
    # the look-up key is a range reference, brs searches an array of range references
    my $locus = binary_range_search(
        [ $gff_line->{start}, $gff_line->{end} ],
        $sorted_annotations{ $gff_line->{seqname} },
    ) || next COORD;

    # Step 4: Reverse search if needed
    # orientation of search: 5'->3' or 5'<-3' (reverse)
    my $reverse
        = ( $five_prime and $locus->[2] eq q{-} )
            || ( $three_prime and $locus->[2] eq q{+} )
                || 0;

    # Step 5: Locate bin of width bin_width where this score fits into
    # the bin index is the locus 5' or 3' coordinate
    # minus the current start coordinate
    # minus bin_width * (num_bins / 2)
    # divided by the negative bin width.
    # note: $reverse + 4 accesses original start or end coordinate for this locus
    my $index = int(
        ($locus->[ $reverse + 4 ] 
         - $bin_width * ( $num_bins / 2 )
         - $gff_line->{start}
        ) / -$bin_width
    );

    # Step 6: Weight score by amount of overlap
    my $score = weight_score( $locus, $gff_line );

    # Step 7: Save score into its appropriate bin, record its strand
    push @{ $genes{ $locus->[3] }->{scores}[$index] }, $score;
    $genes{ $locus->[3] }->{strand} = $locus->[2];
}

# Step 8: Iterate over scores for each gene and print those
my $scores_iterator = make_scores_iterator( \%genes );
while ( my ( $locus, $scores_ref ) = $scores_iterator->($num_bins) ) {
    print join( "\t", $locus, @{$scores_ref} ), "\n";
}


## done


sub weight_score {
    my ( $locus, $gff_line ) = @_;

    my $bin_low  = max( $locus->[0], $gff_line->{start} );
    my $bin_high = min( $locus->[1], $gff_line->{end} );
    my $overlap  = abs($bin_high - $bin_low) + 1;
    my $weight   = $overlap / ( $gff_line->{end} - $gff_line->{start} + 1 );

    return $gff_line->{score} * $weight;
}

sub make_scores_iterator {
    my ($genes_ref) = @_;
    my @genes = sort keys %{$genes_ref};

    return sub {
        my ($num_bins) = @_;

        my $gene = shift @genes || return;
        my @scores = ();

        for my $index ( 0 .. $num_bins - 1 ) {

            # reverse bin orientation for reverse strand genes
            $index = $num_bins - $index - 1
                if $genes{$gene}->{strand} eq q{-};

            push @scores,
                (
                defined $genes{$gene}->{scores}[$index]
                ? sprintf( "%g",
                    ( sum @{ $genes{$gene}->{scores}[$index] } )
                        / @{ $genes{$gene}->{scores}[$index] } )
                : 'na'
                );
        }
        return ( $gene, \@scores );
    }
}

sub binary_range_search {
    my ( $range, $ranges ) = @_;
    my ( $low, $high ) = ( 0, @{$ranges} - 1 );
    while ( $low <= $high ) {

        my $try = int( ( $low + $high ) / 2 );

        $low = $try + 1, next if $ranges->[$try][1] < $range->[0];
        $high = $try - 1, next if $ranges->[$try][0] > $range->[1];

        # range_assert ($range->[0], @{$ranges->[$try]});
        return $ranges->[$try];
    }
    return;
}

sub range_assert {
    my ( $coord, $low, $high ) = @_;

    if ( $coord < $low or $coord > $high ) {
        croak "NOT OK: ", $coord, "\t$low-$high\n";
    }
}

sub make_gff_iterator {
    my ( $gff_file_name, $gff_parser ) = @_;

    open my $GFF_FILE, '<', $gff_file_name
        or croak "Can't read $gff_file_name: $!";

    return sub { $gff_parser->( scalar <$GFF_FILE> ) };
}

sub gff_read {
    return [] if $_[0] =~ m/^
                            \s*
                            \#+
                           /mx;

    my ($seqname, $source, $feature, $start, $end,
        $score,   $strand, $frame,   $attribute
    ) = split m/\t/xm, shift || return;

    $attribute =~ s/[\r\n]//mxg;

    return {
        'seqname'   => lc $seqname,
        'source'    => $source,
        'feature'   => $feature,
        'start'     => $start,
        'end'       => $end,
        'score'     => $score,
        'strand'    => $strand,
        'frame'     => $frame,
        'attribute' => $attribute
    };
}

sub index_gff_annotation {
    my ( $gff_file_name, $attribute_id ) = @_;

    my $gff_annotation = {};
    my $gff_iterator = make_gff_iterator( $gff_file_name, \&gff_read );

 GFF:
    while ( my $gff_line = $gff_iterator->() ) {

        next GFF unless ref $gff_line eq 'HASH';

	$gff_line->{attribute} =~ s/.*
                                    $attribute_id
                                    =
                                    (\w+)
                                    .*
                                   /$1/mx
                                       if $attribute_id;

        $gff_annotation->{ $gff_line->{seqname} }{ $gff_line->{start} } = [
            $gff_line->{start},  $gff_line->{end},
            $gff_line->{strand}, $gff_line->{attribute}
        ];
    }
    return $gff_annotation;
}

sub offset_gff_annotation {
    my ( $gff_annotation, $flag_parser, $parameters ) = @_;
    my $offset_gff_annotation = {};

    return unless $gff_annotation and $flag_parser;

    for my $seqid ( sort keys %{$gff_annotation} ) {

        # every count of 3 genes, we offset the middle one
        # this ensures that we don't miss the first gene
        my @memory = ( [ 0, 0, q{+}, q{} ] );

        for my $start (
            sort { $a <=> $b }
                keys %{ $gff_annotation->{$seqid} }
            ) {

            push @memory, $gff_annotation->{$seqid}{$start};

            # buffer is full, offset middle gene and delete first gene
            if ( @memory == 3 ) {
                check_neighbourhood( $offset_gff_annotation, $flag_parser,
                    $parameters, $seqid, @memory );
                shift @memory;
            }
        }

        # this ensures that we don't miss the last gene
        push @memory, [ 0, 0, q{+}, q{} ];
        check_neighbourhood( $offset_gff_annotation, $flag_parser,
            $parameters, $seqid, @memory );

    }
    return $offset_gff_annotation;
}

sub check_neighbourhood {
    my ( $offset_gff_annotation, $flag_parser, $parameters, $seqid, @memory )
        = @_;

    # offset gene coordinates using appropriate flag method
    my ( $flag_start, $flag_end ) = $flag_parser->( $parameters, @memory );

    $offset_gff_annotation->{$seqid}{$flag_start} = [
        $flag_start,     $flag_end,       $memory[1]->[2],
        $memory[1]->[3], $memory[1]->[0], $memory[1]->[1]
    ]
        if $flag_start <= $flag_end;

    return 1;
}



sub stop_flag_0 {
    my ( $parameters, $previous, $current, $next ) = @_;

    return ( min_max_distance( $current, $parameters ) );
}

sub stop_flag_1 {
    croak "Not yet implemented";
}

sub stop_flag_2 {
    my ( $parameters, $previous, $current, $next ) = @_;

    my ( $minimum, $maximum ) = min_max_distance( $current, $parameters );

    if (   $parameters->{-five_prime} and $current->[2] eq q{-}
        or $parameters->{-three_prime} and $current->[2] eq q{+} )
    {
        return ( max( $current->[0], $previous->[1], $minimum ),
                 min( $next->[0], $maximum ) );
    }
    else {
        return ( max( $previous->[1], $minimum ),
                 min( $current->[1], $next->[0], $maximum ) );
    }
}

sub stop_flag_3 {
    croak "Not yet implemented";
}

sub stop_flag_4 {
    croak "Not yet implemented";
}

sub stop_flag_5 {
    croak "Not yet implemented";
}

sub stop_flag_6 {
    my ( $parameters, $previous, $current, $next ) = @_;

    croak "Flag 6 requires the stop distance parameter (-k)"
        unless $parameters->{-stop_distance};

    my ( $minimum, $maximum ) = min_max_distance( $current, $parameters );

    if (   $parameters->{-five_prime} and $current->[2] eq q{-}
        or $parameters->{-three_prime} and $current->[2] eq q{+} )
    {
        return (
            max($current->[0] + $parameters->{-stop_distance}, $previous->[1],
                $minimum
            ),
            min( $next->[0], $maximum )
        );
    }
    else {
        return (
            max($previous->[1], $minimum ),
            min($current->[1] - $parameters->{-stop_distance}, $next->[0], $maximum)
        );
    }
}

sub min_max_distance {
    my ( $current, $parameters ) = @_;

    my ( $minimum, $maximum );

    if (   $parameters->{-five_prime} and $current->[2] eq q{-}
        or $parameters->{-three_prime} and $current->[2] eq q{+} )
    {

        $minimum = $current->[1] - $parameters->{-distance};
        $maximum = $current->[1] + $parameters->{-distance};

    }
    else {
        $minimum = $current->[0] - $parameters->{-distance};
        $maximum = $current->[0] + $parameters->{-distance};
    }

    $minimum = ( $minimum < 0 ? 0 : $minimum );

    return ( $minimum, $maximum );
}

__END__


=head1 NAME

 ends_analysis.pl - Produce histogram of GFF data scores, given a GFF annotation

=head1 SYNOPSIS

 ends_analysis.pl -g gene_annotation.gff -b 100 -d 5000 -s 2 -5 -x ID probes.gff

=head1 DESCRIPTION

=head1 OPTIONS

 ends_analysis.pl [OPTION]... [FILE]...

 -g, --gff-annotation GFF 3 annotation file
 -b, --bin-width      histogram bin width                                [100]
 -d, --distance       bp distance from end terminal to search, both ways [5000]
 -s, --stop-flag      when to stop searching (FIXME: explain options)    [2]
 -k, --stop-distance  distance from genes to stop from (flag 6 only)     [1500]
 -3, --three-prime    center analysis on 3' end                          [0]
 -5, --five-prime     center analysis on 5' end                          [1]
 -x, --extract-id     GFF 3 attribute tag pointing to locus ID           [ID]
 -o, --output         filename to write results to (defaults to STDOUT)
 -v, --verbose        output perl's diagnostic and warning messages
 -q, --quiet          supress perl's diagnostic and warning messages
 -h, --help           print this information
 -m, --manual         print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev$:
 $Author$:
 $Date$:
 $HeadURL$:
 $Id$:

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
