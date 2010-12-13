#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/min max/;

my $DATA_HANDLE = 'ARGV';
my $gff;
my $tag         = 'ID';
my $full;
my $overlap     = 0.0;
my $no_sort;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'gff|g=s'     => \$gff,
    'tag|t=s'     => \$tag,
    'full|f'      => \$full,
    'no-sort|n'   => \$no_sort,
    'overlap|p=f' => \$overlap,
    'output|o=s'  => \$output,
    'verbose|v'   => sub { use diagnostics; },
    'quiet|q'     => sub { no warnings; },
    'help|h'      => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'    => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

# slurp up GFF's given as arguments.  these will serve as source's.

my $gff_iterator
= make_gff_iterator( parser => \&gff_read, handle => 'ARGV' );

my %gff_records = (); # { sequence => [$gff_hashes] }
LOAD:
while ( my $gff_line = $gff_iterator->() ) {
    next LOAD unless ref $gff_line eq 'HASH';
    push @{ $gff_records{ $gff_line->{seqname} } }, $gff_line;
}

# $gff file will serve as target's

my $annotation_iterator = sub { };
if ($gff) {
    $annotation_iterator = make_annotation_iterator(
        file  => $gff,
        tag   => $tag,
    );
}

SEQUENCE:
for my $sequence ( sort keys %gff_records ) {

    # sort gff records by start (maybe)
    unless ($no_sort) {
        @{ $gff_records{$sequence} }
            = sort { $a->{start} <=> $b->{start} }
            @{ $gff_records{$sequence} };
    }

    # for each range/locus for the sequence...
  WINDOW:
    while ( my ( $ranges, $locus ) = $annotation_iterator->($sequence) ) {


        my $brs_iterator = binary_range_search(
            # $ranges specific for a particular $locus
            range  => $ranges,
            # all $ranges for the sequence, sorted by start
            ranges => $gff_records{$sequence},
        );

        # what is $_ here???
        next WINDOW unless
        $overlap <= overlap ($ranges->[0], [$_->{start}, $_->{end}]);

        # print if overlap meets a cutoff.
        if ($full) {
            $overlap <= overlap ($ranges->[0], [$_->{start}, $_->{end}])
            && print join ("\t",
                           $_->{seqname},
                           $_->{source},
                           $_->{feature},
                           $_->{start},
                           $_->{end},
                           $_->{score},
                           $_->{strand},
                           $_->{frame},
                           $_->{attribute},
                       ), "\n"
                       while local $_ = $brs_iterator->();
        }
        else {
            my $overlaps = 0;

            $overlap <= overlap ($ranges->[0], [$_->{start}, $_->{end}])
            && $overlaps++ while local $_ = $brs_iterator->();

            print join ("\t", 
                        $locus,
                        $overlaps
                    ), "\n";
        }
    }
}

##########################################################
# overlap
#
# a: |--------------|
# b:    |----------------------------|
#       |-----------| <- overlap region
# 
# return ratio of overlap region length / range_b length

sub overlap {
    
    my ($range_a, $range_b) = @_;

    # the higher start range
    my $bin_low  = max( $range_a->[0], $range_b->[0] );

    # the lower end range 
    my $bin_high = min( $range_a->[1],   $range_b->[1]   );

    # length of overlap
    my $overlap  = abs($bin_high - $bin_low) + 1;

    # overlap ratio = overlap / length of range_b;
    return $overlap / ( $range_b->[1] - $range_b->[0] + 1 );
}

##########################################################
# make_annotation_iterator
#
# return an iterator which, given a sequence, will return 
# ( [ ranges ], $locus ), on every call, until loci runs out

sub make_annotation_iterator {
    my (%options) = @_;

    my $annotation_file = $options{file};
    my $locus_tag       = $options{tag}   || 'ID';
    my $merge_exons     = $options{merge} || 0;

    open my $GFFH, '<', $annotation_file
        or croak "Can't read file: $annotation_file";

    my $gff_iterator
        = make_gff_iterator( parser => \&gff_read, handle => $GFFH );

    my $annotation = index_annotation( iterator => $gff_iterator, %options );

    close $GFFH or croak "Can't close $annotation_file: $!";

    # closure
    my %annotation_keys = ();
    return sub {
        my ($sequence) = @_;
        return unless $sequence;

        # first time, make %annotation_keys a { sequence => [locus_ids] }
        if ( ! exists $annotation_keys{$sequence} ) {
            @{ $annotation_keys{$sequence} }
                = sort keys %{ $annotation->{$sequence} };
        }

        # grab a locus
        my $locus = shift @{ $annotation_keys{$sequence} };
        return unless $locus;

        # grab the locus's [ranges]
        my $ranges = $annotation->{$sequence}{$locus};
        delete $annotation->{$sequence}{$locus};

        # return [ unique_ranges of seq ], and a locus
        if ( $locus and @$ranges ) {
            return [ uniq_ranges($ranges) ], $locus;
        }
        else {
            return;
        }
    };
}

##########################################################
# index annotations
# 
# accept a gff iterator, a locus tag ('ID') and a feature (column 3) to merge
# on.
#
# return { sequence_names => { locus_id => [ranges]} }

sub index_annotation {
    my (%options) = @_;

    my $gff_iterator = $options{iterator}
        || croak 'Need an iterator parameter';
    my $locus_tag     = $options{tag}   || 'ID';
    my $merge_feature = $options{merge} || 0;

    my %annotation = ();

LOCUS:
    # fore each line in gff
    while ( my $locus = $gff_iterator->() ) {

        next LOCUS unless ref $locus eq 'HASH';

        # see if there's an 'ID'. 
        my ($locus_id) = $locus->{attribute} =~ m/$locus_tag[=\s]?([^;,]+)/;

        # if there is no ID, use the first attribute 
        if ( !defined $locus_id ) {
            ( $locus_id, undef ) = split /;/, $locus->{attribute};
            $locus_id ||= q{.};
        }
        # if there is...
        else {

            # clean whitespace/quotes
            $locus_id =~ s/["\t\r\n]//g;

            # kill after period if merging
            $locus_id
                =~ s/\.\w+$// # this fetches the parent ID in GFF gene models (eg. exon is Parent=ATG101010.n)
                if $merge_feature and $locus->{feature} eq $merge_feature;
        }

        # so annotation is a { sequence_names => { locus_id => [ranges]} }
        push @{ $annotation{ $locus->{seqname} }{$locus_id} },
            [ $locus->{start}, $locus->{end} ]
            unless ( $merge_feature and $locus->{feature} ne $merge_feature );

    }

    return \%annotation;
}

##########################################################
# remove duplicate ranges
# 
# $ranges is an array ref of ranges
# [[$start_a, $end_a], [$start_b, $end_b], ... ]

sub uniq_ranges {
    my ($ranges) = @_;

    my %seen;
    my @uniq;

    for (@$ranges) {
        push @uniq, $_ unless $seen{"$_->[0];$_->[1]"}++;
    }

    return wantarray ? @uniq : [@uniq];
}

##########################################################
# binary_range_search
# 
# given a list of target ranges, and a list of source ranges,
# return an iterator which iterates over source ranges that overlaps with
# a target range
#

sub binary_range_search {
    my %options = @_;

    # list of ranges for a particular locus and sequence
    my $targets = $options{range}  || croak 'Need a range parameter';

    # sorted (by range-starts) of all ranges for sequence
    my $ranges  = $options{ranges} || croak 'Need a ranges parameter';

    # lowest and highest indicies of $ranges
    my ( $low, $high ) = ( 0, $#{$ranges} );
    my @iterators = ();

TARGET:
    for my $range (@$targets) {

    RANGE_CHECK:
        while ( $low <= $high ) {

            # middle
            my $try = int( ( $low + $high ) / 2 );

            # try higher $try if:
            # $range: |------------| 
            # $try:                  |---------------|
            $low = $try + 1, next RANGE_CHECK
                if $ranges->[$try]{end} < $range->[0];

            # try lower $try if:
            # $range:                   |------------| 
            # $try:   |---------------|
            # try lower half if range's end < middle's start
            $high = $try - 1, next RANGE_CHECK
                if $ranges->[$try]{start} > $range->[1];

            # if we make it here, we've found a source $try which overlaps with the  
            # target $range
            
            my ( $down, $up ) = ($try) x 2;
            my %seen = ();

            # create an iterator which, on every call, returns an overlapping
            # source range, starting with $range and then going down or up the
            # list.
            my $brs_iterator = sub {
                
                # if $range overlaps with $up + 1, and it's new, return it
                if (    $ranges->[ $up + 1 ]{end} >= $range->[0]
                    and $ranges->[ $up + 1 ]{start} <= $range->[1]
                    and !exists $seen{ $up + 1 } )
                {
                    $seen{ $up + 1 } = undef;
                    return $ranges->[ ++$up ];
                }
                # if $range overlaps with $down - 1, and it's new, return it
                elsif ( $ranges->[ $down - 1 ]{end} >= $range->[0]
                    and $ranges->[ $down - 1 ]{start} <= $range->[1]
                    and !exists $seen{ $down - 1 }
                    and $down > 0 )
                {
                    $seen{ $down - 1 } = undef;
                    return $ranges->[ --$down ];
                }
                # we already know $try overlaps, so return it too.
                elsif ( !exists $seen{$try} ) {
                    $seen{$try} = undef;
                    return $ranges->[$try];
                }
                else {
                    return;
                }
            };
            push @iterators, $brs_iterator;
            next TARGET;
        }
    }

    # In scalar context return master iterator that iterates over the list of range iterators.
    # In list context returns a list of range iterators.
    return wantarray
        ? @iterators
        : sub {
        while (@iterators) {
            if ( my $range = $iterators[0]->() ) {
                return $range;
            }
            shift @iterators;
        }
        return;
        };
}

##########################################################
# Gff reading, to be moved out

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

sub make_gff_iterator {
    my %options = @_;

    my $parser     = $options{parser};
    my $file       = $options{file};
    my $GFF_HANDLE = $options{handle};

    croak
        "Need parser function reference and file name or handle to build iterator"
        unless $parser
            and ref $parser eq 'CODE'
            and (
                ( defined $file and -e $file )
                xor(defined $GFF_HANDLE and ref $GFF_HANDLE eq 'GLOB'
                        or $GFF_HANDLE eq 'ARGV'
                )
            );

    if ($file) {
        open $GFF_HANDLE, '<', $file
            or croak "Can't read $file: $!";
    }

    return sub {
        $parser->( scalar <$GFF_HANDLE> );
    };
}


__END__


=head1 NAME

 overlaps_gff.pl - Compute overlaps on GFF data and annotations

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

 overlaps_gff.pl [OPTION]... [FILE]...

 -g, --gff         GFFv3 annotation file on which to look for overlaps
 -t, --tag         GFFv3 attribute tag ('ID') on which value to look for loci IDs
 -f, --full        Don't just print loci and overlap number
 -p, --overlap     accepted overlapping coefficient
 -n, --no-sort     Input GFFv3 data is sorted numerically
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
