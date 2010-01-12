#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
Getopt::Long::Configure('bundling');
use Pod::Usage;

my $output;
my $width   = 50;
my $step    = 50;
my $merge   = 0;
my $no_sort = 0;
my $no_skip = 0;
my $average = 0;
my $gff;
my $tag     = 'ID';

# Grabs and parses command line options
my $result = GetOptions(
    'width|w=i'  => \$width,
    'step|s=i'   => \$step,
    'average|a'  => \$average,
    'merge|m'    => \$merge,
    'no-sort|n'  => \$no_sort,
    'no-skip|k'  => \$no_skip,
    'gff|g=s'    => \$gff,
    'tag|t=s'    => \$tag,
    'output|o=s' => \$output,
    'verbose'    => sub { use diagnostics; },
    'quiet'      => sub { no warnings; },
    'help'       => sub { pod2usage( -verbose => 1 ); },
    'manual'     => sub { pod2usage( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage( -verbose => 1 )
    unless @ARGV and $result;

if ($output) {
    open my $USER_OUT, '>', $output
        or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

if ($merge) {
    $width = 1;
    $step  = 1;
}

my $gff_iterator
    = make_gff_iterator( parser => \&gff_read, handle => 'ARGV' );

my %gff_records = ();

LOAD:
while ( my $gff_line = $gff_iterator->() ) {

    # &gff_read returns [] for GFF comments, invalid lines, etc.
    next LOAD unless ref $gff_line eq 'HASH';

    push @{ $gff_records{ $gff_line->{seqname} } }, $gff_line;

}

my $window_iterator;

if ($gff) {
    $window_iterator = make_window_iterator(
        gff => {
            file  => $gff, 
            tag   => $tag, 
            merge => $merge
        }
    );
}

SEQUENCE:
for my $sequence ( sort keys %gff_records ) {

    unless ($no_sort) {
        @{ $gff_records{$sequence} } = sort { $a->{start} <=> $b->{start} }
            @{ $gff_records{$sequence} };
    }

    unless ($gff) {
        $window_iterator = make_window_iterator(
            windows => {
                width => $width,
                step  => $step, 
                lower => 1,
                upper => $gff_records{$sequence}[-1]->{end}
            }
        );
    }

  WINDOW:
    while ( my ($ranges, $locus) = $window_iterator->($sequence) ) {

        # print STDERR "$start\t$end\n";next WINDOW;

        my $brs_iterator = binary_range_search(
            range    => $ranges,
            ranges   => $gff_records{$sequence},
        );

        my $scores_ref
            = $average
            ? average_scores($brs_iterator)
            : fractional_methylation($brs_iterator);

        my ( $score, $attribute ) = qw(. .);

        if ( ref $scores_ref eq 'HASH' ) {

            $score = sprintf( "%g", $scores_ref->{score} );
            delete $scores_ref->{score};

            $attribute = join q{; }, map { "$_=" . $scores_ref->{$_} }
                sort keys %{$scores_ref};

        }

        print join( "\t",
                    $sequence, 'dzlab',         "w$width",
                    $start,        $end, $score,
                    q{.},      q{.},            $attribute,
                ),
                "\n";

    }

    delete $gff_records{$sequence};
}

sub make_window_iterator {
    my (%options) = @_;

    if ($options{windows}) {
        my $width = $options{windows}{width};
        my $step  = $options{windows}{step};
        my $lower = $options{windows}{lower};
        my $upper = $options{windows}{upper};
        
        return sub {
            my $i      = $lower;
            $lower    += $step;

            if ($i <= $upper - $width + 1) {
                return [ [$i, $i + $width - 1] ], "w$step";
            }
            else {
                return undef;
            }
        }
    }
    elsif ($options{gff}) {

        my $annotation_file = $options{gff}{file};
        my $locus_tag       = $options{gff}{tag}   || 'ID';
        my $merge_exons     = $options{gff}{merge} || 0;
        
        open my $GFFH, '<', $annotation_file 
        or croak "Can't read file: $annotation_file";

        my $gff_iterator
        = make_gff_iterator( parser => \&gff_read, handle => $GFFH );

        my %annotation = ();

      LOCUS:
        while ( my $locus = $gff_iterator->() ) {

            next LOCUS unless ref $locus eq 'HASH';

            my ($locus_id) = $locus->{attribute} =~ m/$locus_tag[=\s]?([^;]+)/;

            if (!defined $locus_id) {
                ($locus_id, undef) = split /;/, $locus->{attribute};
                $locus_id ||= q{.};
            } 
            else {
                $locus_id =~ s/["\t\r\n]//g;
                $locus_id =~ s/\.\d$// if $merge_exons;
            }

            # push @{ $annotation{$locus->{seqname}}{$locus->{start}}},
            # [$locus->{start}, $locus->{end}, $locus_id, $locus->{strand}, $locus->{source}, $locus->{feature}, $locus->{attribute}];

            push @{ $annotation{$locus->{seqname}}{$locus_id} },
            [$locus->{start}, $locus->{end}];

        }
        close $GFFH;

        my %annotation_keys;

        return sub {
            my ($sequence) = @_;
            return unless $sequence;

            unless (exists $annotation_keys{$sequence}) {
                @{ $annotation_keys{$sequence} }
                = sort {$annotation{$sequence}{$a}->[0] <=> $annotation{$sequence}{$b}->[0]}
                    keys %{$annotation{$sequence}};
            }

            my $locus  = shift @{ $annotation_keys{$sequence} };
            my $ranges = shift @{ $annotation{$sequence}{$locus} };

            if ($locus and @$ranges) {
                return $ranges, $locus;
            }
            else {
                return undef;
            }
        };
    }
}


sub fractional_methylation {
    my ($brs_iterator) = @_;

    my ( $c_count, $t_count, $score_count ) = ( 0, 0, 0 );

COORD:
    while ( my $gff_line = $brs_iterator->() ) {

        next COORD unless ref $gff_line eq 'HASH';

        my ( $c, $t ) = $gff_line->{attribute} =~ m/
                                                     c=(\d+)
                                                     .*
                                                     t=(\d+)
                                                 /xms;

        next unless defined $c and defined $t;

        $c_count += $c;
        $t_count += $t;
        $score_count++;
    }

    if ( $score_count and ( $c_count or $t_count ) ) {
        return {
            score => $c_count / ( $c_count + $t_count ),
            c     => $c_count,
            t     => $t_count,
            n     => $score_count,
        };
    }
}

sub average_scores {
    my ($brs_iterator) = @_;

    my ( $score_avg, $score_std, $score_var, $score_count ) = ( 0, 0, 0, 0 );

COORD:
    while ( my $gff_line = $brs_iterator->() ) {

        next COORD unless ref $gff_line eq 'HASH';

        my $previous_score_avg = $score_avg;

        $score_avg = $score_avg
            + ( $gff_line->{score} - $score_avg ) / ++$score_count;

        $score_std
            = $score_std
            + ( $gff_line->{score} - $previous_score_avg )
            * ( $gff_line->{score} - $score_avg );

        $score_var = $score_std / ( $score_count - 1 )
            if $score_count > 1;

    }

    if ($score_count) {
        return {
            score => $score_avg,
            std   => sqrt ($score_var),
            var   => $score_var,
            n     => $score_count,
        };
    }
}

sub range_assert {
    my ( $coord, $low, $high ) = @_;

    if ( $coord < $low or $coord > $high ) {
        croak "NOT OK: $coord\t$low-$high\n";
    }
}

sub binary_range_search {
    my %options = @_;
    my $rage     = $options{range}  || return;
    my $ranges   = $options{ranges} || return;

    my ( $low, $high ) = ( 0, @{$ranges} - 1 );

    while ( $low <= $high ) {

        my $try = int( ( $low + $high ) / 2 );

        $low  = $try + 1, next if $ranges->[$try]{end}   < $range->[0];
        $high = $try - 1, next if $ranges->[$try]{start} > $range->[1];

        my ( $down, $up ) = ($try) x 2;

        my %seen = ();

        my $brs_iterator = sub {

            if (    $ranges->[ $up + 1 ]{end}       >= $range->[0]
                    and $ranges->[ $up + 1 ]{start} <= $range->[1]
                    and !exists $seen{ $up + 1 } )
            {
                $seen{ $up + 1 } = undef;
                return $ranges->[ ++$up ];
            }
            elsif ( $ranges->[ $down - 1 ]{end}       >= $range->[0]
                    and $ranges->[ $down + 1 ]{start} <= $range->[1]
                    and !exists $seen{ $down - 1 }
                    and $down > 0 )
            {
                $seen{ $down - 1 } = undef;
                return $ranges->[ --$down ];
            }
            elsif ( !exists $seen{$try} ) {
                $seen{$try} = undef;
                return $ranges->[$try];
            }
            else {
                return;
            }

        };
        return $brs_iterator;
    }
    return sub { };
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

 $Rev$:
 $Author$:
 $Date$:
 $HeadURL$:
 $Id$:

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
222:	final indentation level: 1

Final nesting depth of '{'s is 1
The most recent un-matched '{' is on line 87
87: sub binary_range_search {
                            ^
222:	To save a full .LOG file rerun with -g
