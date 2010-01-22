#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
Getopt::Long::Configure('bundling');
use Pod::Usage;

my $output;
my $width    = 50;
my $step     = 50;
my $merge    = 0;
my $no_sort  = 0;
my $no_skip  = 0;
my $scoring  = 'meth'; # meth, average, or sum
my $gff;
my $tag      = 'ID';
my $feature;
my $absolute = 0;

# Grabs and parses command line options
my $result = GetOptions(
    'width|w=i'    => \$width,
    'step|s=i'     => \$step,
    'scoring|c=s'  => \$scoring,
    'merge|m=s'    => \$merge,
    'no-sort|n'    => \$no_sort,
    'no-skip|k'    => \$no_skip,
    'gff|g=s'      => \$gff,
    'tag|t=s'      => \$tag,
    'feature|f=s'  => \$feature,
    'absolute|b=s' => \$absolute,
    'output|o=s'   => \$output,
    'verbose'      => sub { use diagnostics; },
    'quiet'        => sub { no warnings; },
    'help'         => sub { pod2usage( -verbose => 1 ); },
    'manual'       => sub { pod2usage( -verbose => 2 ); }
);

my %scoring_dispatch = (
    meth    => \&fractional_methylation,
    average => \&average_scores,
    sum     => \&sum_scores,
);

# Check required command line parameters
pod2usage( -verbose => 1 )
unless @ARGV and $result
and (($width and $step) or $gff)
and exists $scoring_dispatch{$scoring};

if ($output) {
    open my $USER_OUT, '>', $output
        or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

if ($merge) {
    $width = 1;
    $step  = 1;
}

if ($absolute) {
    croak '-b, --absolute option only works with arabidopsis or rice'
    unless lc $absolute eq 'arabidopsis' or lc $absolute eq 'rice' or lc $absolute eq 'puffer';
    $no_skip = 1;
}

my $gff_iterator
    = make_gff_iterator( parser => \&gff_read, handle => 'ARGV' );

my %gff_records = ();
LOAD:
while ( my $gff_line = $gff_iterator->() ) {
    next LOAD unless ref $gff_line eq 'HASH';
    push @{ $gff_records{ $gff_line->{seqname} } }, $gff_line;
}

my $fasta_lengths = {};
if ($absolute and $no_skip) {
    $fasta_lengths
    = index_fasta_lengths( handle => 'DATA' );
}

my $window_iterator = sub {};
if ($gff) {
    $window_iterator = make_annotation_iterator(
        file  => $gff, 
        tag   => $tag, 
        merge => $merge
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
            width => $width,
            step  => $step, 
            lower => 1,
            upper => ($absolute and $no_skip and %$fasta_lengths
                      and exists $fasta_lengths->{$absolute}{$sequence})
                      ? $fasta_lengths->{$absolute}{$sequence}
                      : $gff_records{$sequence}[-1]->{end},
        );
    }

  WINDOW:
    while ( my ($ranges, $locus) = $window_iterator->($sequence) ) {

        my $brs_iterator = binary_range_search(
            range    => [@$ranges], # make an anonymous copy because binary_range_search *will* modify references
            ranges   => $gff_records{$sequence},
        );

        my $scores_ref
        = $scoring_dispatch{$scoring}->($brs_iterator);

        my ( $score, $attribute );

        if ( ref $scores_ref eq 'HASH' ) {

            $score = sprintf( "%g", $scores_ref->{score} );
            delete $scores_ref->{score};

            $attribute  = "ID=$locus; " if $locus;
            $attribute .= join q{; }, map { "$_=" . sprintf("%g", $scores_ref->{$_}) }
                sort keys %{$scores_ref};

        }
        elsif ($no_skip) {
            $score     = q{.};
            $attribute = $locus ? "ID=$locus" : q{.};
        }
        else {
            next WINDOW;
        }

        print join( "\t",
                    $sequence, 'dzlab', ($feature ? $feature : $gff ? 'locus' : 'window'),
                    $ranges->[0][0], $ranges->[-1][1], $score,
                    q{.}, q{.}, $attribute,
                ), "\n";
    }

    delete $gff_records{$sequence};
}


sub index_fasta_lengths {
    my %options = @_;

    my $handle = $options{handle};

    if ($options{file}) {
        open $handle, '<', $options{file} or croak $!;
    }

    my %fasta_lengths;
    my  $active;
    while (<$handle>) {
        chomp;
        if (m/^\s*#/) {
            s/#//;
            $active = lc $_;
        }
        else {
            my ($sequence, $length) = split /\t/;
            $fasta_lengths{$active}{$sequence} = $length;
        }
    }

    if ($options{file}) {
        close $handle or croak $!;
    }

    return \%fasta_lengths;
}


sub make_window_iterator {
    my (%options) = @_;

    my $width = $options{width} || croak 'Need window parameter';
    my $step  = $options{step}  || croak 'Need step parameter';
    my $lower = $options{lower} || croak 'Need lower bound parameter';
    my $upper = $options{upper} || croak 'Need upper bound parameter';
        
    return sub {
        my $i      = $lower;
        $lower    += $step;

        if ($i <= $upper - $width + 1) {
            return [ [$i, $i + $width - 1] ];
        } 
        elsif ($i < $upper) {
            return [ [$i, $upper] ];
        }
        else {
            return;
        }
    }
}

sub make_annotation_iterator {
    my (%options) = @_;

    my $annotation_file = $options{file};
    my $locus_tag       = $options{tag}   || 'ID';
    my $merge_exons     = $options{merge} || 0;
        
    open my $GFFH, '<', $annotation_file 
    or croak "Can't read file: $annotation_file";

    my $gff_iterator
    = make_gff_iterator( parser => \&gff_read, handle => $GFFH );

    my $annotation 
    = index_annotation( iterator => $gff_iterator, %options );

    close $GFFH or croak "Can't close $annotation_file: $!";

    my %annotation_keys = ();
    return sub {
        my ($sequence) = @_;
        return unless $sequence;

        unless (exists $annotation_keys{$sequence}) {
            @{ $annotation_keys{$sequence} }
            = sort keys %{$annotation->{$sequence}};
        }

        my $locus  = shift @{ $annotation_keys{$sequence} };
        return unless $locus;
        my $ranges = $annotation->{$sequence}{$locus}; delete $annotation->{$sequence}{$locus};

        if ($locus and @$ranges) {
            return [uniq_ranges ($ranges)], $locus;
        } 
        else {
            return;
        }
    };
}


sub index_annotation {
    my (%options) = @_;

    my $gff_iterator  = $options{iterator} || croak 'Need an iterator parameter';
    my $locus_tag     = $options{tag}      || 'ID';
    my $merge_feature = $options{merge}    || 0;

    my %annotation    = ();

  LOCUS:
    while ( my $locus = $gff_iterator->() ) {

        next LOCUS unless ref $locus eq 'HASH';

        my ($locus_id) = $locus->{attribute} =~ m/$locus_tag[=\s]?([^;,]+)/;

        if (!defined $locus_id) {
            ($locus_id, undef) = split /;/, $locus->{attribute};
            $locus_id ||= q{.};
        } 
        else {
            $locus_id =~ s/["\t\r\n]//g;
            $locus_id =~ s/\.\d$// # this fetches the parent ID in GFF gene models (eg. exon is Parent=ATG101010.n)
            if $merge_feature and $locus->{feature} eq $merge_feature;
        }

        push @{ $annotation{$locus->{seqname}}{$locus_id} },
        [$locus->{start}, $locus->{end}] 
        unless ($merge_feature and $locus->{feature} ne $merge_feature);

    }

    return \%annotation;
}



sub uniq_ranges {
    my ($ranges) = @_;

    my %seen;
    my @uniq;

    for (@$ranges) {
        push @uniq, $_ unless $seen{"$_->[0];$_->[1]"}++;
    }

    return wantarray ? @uniq : [@uniq];
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

sub sum_scores {
    my ($brs_iterator) = @_;

    my ($score_sum, $score_count) = (0, 0);

  COORD:
    while ( my $gff_line = $brs_iterator->() ) {
        next COORD unless ref $gff_line eq 'HASH';

        $score_sum  += $gff_line->{score} eq q{.} ? 1 : $gff_line->{score};
        $score_count++;
    }

    if ($score_count) {
        return {
            score => $score_sum,
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


sub binary_range_search {
    my %options = @_;

    my $range  = shift @{$options{range}} || croak 'Need a range parameter';
    my $ranges = $options{ranges}         || croak 'Need a ranges parameter';

    my ( $low, $high ) = ( 0, @{$ranges} - 1 );

    while ( $low <= $high ) {

        my $try = int( ( $low + $high ) / 2 );

        $low  = $try + 1, next if $ranges->[$try]{end}   < $range->[0];
        $high = $try - 1, next if $ranges->[$try]{start} > $range->[1];

        my ( $down, $up ) = ($try) x 2;

        my %seen = ();

        my $brs_iterator = sub { };
        $brs_iterator = sub {

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
            elsif (@{$options{range}}) {
                $brs_iterator = 
                binary_range_search(
                    range    => $options{range},
                    ranges   => $options{ranges},
                );
                return $brs_iterator->();
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
        and ( ( defined $file and -e $file )
              xor (defined $GFF_HANDLE and ref $GFF_HANDLE eq 'GLOB'
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

=head1 NAME

 window_gff.pl - Average GFFv3 data over a sliding window or against a GFFv3 annotation file

=head1 SYNOPSIS

 # Run a sliding window of 50bp every 25bp interval on methylation data (attribute field assumed to be 'c=n; t=m')
 window_gff.pl --width 50 --step 25 --scoring meth --output out.gff in.gff

 # Average data per each locus (on score field)
 window_gff.pl --gff genes.gff --scoring average --output out.gff in.gff

 # Merge 'exon' features in GFF annotation file per parent ID (eg. per gene)
 window_gff.pl --gff exons.gff --scoring average --output out.gff --tag Parent --merge exon

=head1 DESCRIPTION

 This program works in two modes. It will either run a sliding window through the input GFF data,
 or load a GFF annotation file and use that to average scores in input GFF data.

 Scores calculation by default is done via fractional methylation, in which case the program
 requires that the GFF attribute field contain both 'c' and 't' tags (eg. 'c=n; t=m').
 Alternatively, the -a, --average switch will make the program average the scores in the score field.

 An additional switch, -m, --merge, will take the name of a feature in the annotation GFF file.
 It will then try to average out all the features that have a common parent locus.
 This assumes that, for example if averaging exon features: the user must use '--merge exon',
 '--tag Parent' (in the case of Arabidopsis' annotations), which will fetch a gene ID like ATG01010.1.
 It is assumed that the Parent feature of this exon is ATG010101.

=head1 OPTIONS

 window_gff.pl [OPTION]... [FILE]...
 
 -w, --width       sliding window width                                  (default: 50, integer)
 -s, --step        sliding window interval                               (default: 50, integer)
 -c, --scoring     score computation scheme                              (default: meth, string [available: meth, average, sum])
 -m, --merge       merge this feature as belonging to same locus         (default: no, string [eg: exon])
 -n, --no-sort     GFFv3 data assumed sorted by start coordinate         (default: no)
 -k, --no-skip     print windows or loci for which there is no coverage  (deftaul: no)
 -g, --gff         GFFv3 annotation file
 -t, --tag         attribute field tag from which to extract locus ID    (default: ID, string)
 -f, --feature     overwrite GFF feature field with this label           (default: no, string)
 -b, --absolute    organism name to fetch chromosome lengths, implies -k (default: no, string [available: arabidopsis, rice, puffer])
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


__DATA__
#arabidopsis
chr1	30432563
chr2	19705359
chr3	23470805
chr4	18585042
chr5	26992728
chrc	154478
chrm	366924
#rice
chr01	43596771
chr02	35925388
chr03	36345490
chr04	35244269
chr05	29874162
chr06	31246789
chr07	29688601
chr08	28309179
chr09	23011239
chr10	22876596
chr11	28462103
chr12	27497214
chrc	134525
chrm	490520
#puffer
1	22981688
10	13272281
11	11954808
12	12622881
13	13302670
14	10246949
15	7320470
15_random	3234215
16	9031048
17	12136232
18	11077504
19	7272499
1_random	1180980
2	21591555
20	3798727
21	5834722
21_random	3311323
2_random	2082209
3	15489435
4	9874776
5	13390619
6	7024381
7	11693588
8	10512681
9	10554956
mt	16462
