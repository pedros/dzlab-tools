#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;disable diagnostics;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Pod::Usage;

# Globals, passed as command line options
my $ratio_file = q{-}; # gff file with expression, binding sites, etc. data
my $annotation_file;   # gff annotation file
my $distance = 0;      # distance from center of each probe/window on each side to search
my $overlap  = 0;
my $distance_overlap = 0;
my $output   = q{-};
my $verbose  = 0;
my @filter;
my $exclude_non_genes = 0;

# Grabs and parses command line options
my $result = GetOptions (
    'ratio-file|f=s'        => \$ratio_file,
    'annotation-file|a=s'   => \$annotation_file,
    'distance|d=i'          => \$distance,
    'overlap|p'             => \$overlap,
    'distance-overlap|dp|i' => \$distance_overlap,
    'filter|w=i{1,2}'       => \@filter,
    'exclude-non-genes|e'   => \$exclude_non_genes,
    'output|o:s'            => \$output,
    'verbose|v'             => sub {enable diagnostics;use warnings;},
    'quiet|q'               => sub {disable diagnostics;no warnings;},
    'help|usage|h'          => sub {pod2usage(-verbose => 1);},
    'manual|man|m'          => sub {pod2usage(-verbose => 2);}
);

# Check required command line parameters
pod2usage(-verbose => 1)
unless $result && $ratio_file && $annotation_file and ($distance xor $overlap xor $distance_overlap);

# redirects STDOUT to file if specified by user
open STDOUT, '>', "$output" or croak "Can't redirect STDOUT to file: $output"
if $output ne q{-};

# index the annotation file into hash of arrays with hash keys = sequence ids (chromosomes)
# and each array with references to all loci in sequence id in the form [start, end, id, strand]
my %annotation = %{ index_gff_annotation ($annotation_file) };

# slurp expression/binding site data, cleaning up comments, blank lines and de-capitalizing in the process
my %data = (); # %data is a hash of arrays of scalars. Each hash key is a sequence id (chromosome) array reference to a gff line
open my $RATIO, '<', $ratio_file or croak "Can't read file: $ratio_file"
if $ratio_file ne q{-};

while (<$RATIO>) {
    chomp; # delete line feeds
    next if ($_ =~ m/^#.*$|^\s*$|\.;\.$/);
    my $chr = (split /\t/, $_)[0];
    $chr =~ tr/A-Z/a-z/;
    push @{$data{$chr}}, $_;
}
close $RATIO or croak "Can't close file: $ratio_file";

my %logic_dispatch = (
    distance         => sub {
        my $center = int (($_[1] - $_[0]) / 2 + $_[0]);
        return ($center - int $_[2], $center + int $_[2]);
    },
    
    overlap          => sub {
        return ($_[0], $_[1]);
    },

    distance_overlap => sub {
        return ($_[0] - int $_[2], $_[1] + int $_[2]);
    },
);
        
# main logic: process each (sorted) sequence id in turn
# for each one, process each (sorted) window/probe
# expects to find gff formatted line output from chipotle and pre-processed
# to re-arrange columns into standard gff format (chi2gff.pl)
for my $chr (sort {$a cmp $b} keys %data) {

    for my $window (sort {(split /\t/, $a)[3] <=> (split /\t/, $b)[3]} @{$data{$chr}}) {

        my ($feature, $start, $end, $mean, $attribute)
        = (split /\t/, $window)[2, 3, 4, 5, 8];

        if (@filter) {
            my $length = $end - $start + 1;

            next unless
            $length >= $filter[0]
            and ( ($filter[1] and $length <= $filter[1])
                  or not $filter[1]);
        }

        my $center = int (($end - $start) / 2 + $start);
        my $lower_bound;
        my $upper_bound;

        # assume normal distribution within a window (ie peaks are in center of window)
        if ($distance) {
            ($lower_bound, $upper_bound)
            = $logic_dispatch{distance}->($start, $end, $distance);
        }
        # don't assume anything, use overlaps to find genes
        elsif ($overlap) {
            ($lower_bound, $upper_bound)
            = $logic_dispatch{overlap}->($start, $end);
        }
        # go a little $distance beyond the window bounds
        elsif ($distance_overlap) {
            ($lower_bound, $upper_bound)
            = $logic_dispatch{distance_overlap}->($start, $end, $distance_overlap);
        }

        # grab all the loci from the annotation file that fall within that region
        # this is a first filtering step, very coarse. see gff_filter_by_coord
        my @range
        = @{ gff_filter_by_coord (
            $lower_bound, $upper_bound,
            $annotation{$chr}, $overlap
        ) };

        next if
        $exclude_non_genes and not @range;

        # for current window/probe, go through each locus and
        # arrange its parameters (id, distance to center of probe, direction)
        # to put in gff attribute field (form: gene=GENEID:DISTANCE:DIRECTION)
        my $gene_list = @range ? q{gene=} : q{.}; # if no targets, just output '.'
        for my $locus (@range) {

            $gene_list .= $locus->[2]; # start with id

            if ($locus->[0] <= $center && $locus->[1] >= $center) {
                $gene_list .= ':-'; # add '-' if locus overlaps probe center
            }
            else {$gene_list .= ':'}

            if ($locus->[3] eq q{+}) { # check direction of locus
                my $p5 = abs ($locus->[0] - $center);
                my $p3 = abs ($locus->[1] - $center);
                $gene_list .= $p5 <= $p3 ? "$p5:5p" : "$p3:3p";
            }
            else {
                my $p3 = abs ($locus->[0] - $center);
                my $p5 = abs ($locus->[1] - $center);
                $gene_list .= $p5 <= $p3 ? "$p5:5p" : "$p3:3p";
            }
            $gene_list .= q{,}; # gene list separator
        }
        $gene_list =~ s/,$// if @range; # delete the last comma

        # final post-processing to resolve multiple targets
        # logic: if distances between loci are small (<100)
        # go with the 5 prime end, otherwise just go with smaller distance
        # if both candidates are 5prime, go with distance
        if (@range) {
            my @genes
            = split /,/, (split /=/, $gene_list)[1];

            $genes[$_] = [ split /:/, $genes[$_] ]
            for (0 .. @genes - 1);

            my ($min_id, $min_dist, $min_dir)
            = @{$genes[0]}; # initialize best target with first gene

            for my $idx (1 .. @genes - 1) {

                # go with minimal distance
                if ( abs ($genes[$idx]->[1] - $min_dist) > 100) {
                    ($min_id, $min_dist, $min_dir)
                    = @{$genes[$idx]}
                    if $genes[$idx]->[1] < $min_dist;
                }
                else { # or with direction
                    ($min_id, $min_dist, $min_dir)
                    = @{$genes[$idx]}
                    if $genes[$idx]->[2] eq '5p'
                    && ($min_dir ne '5p' or $genes[$idx]->[1] < $min_dist);
                }
            }

            # reconstruct the gene list gff attribute with '*' for best matches
            $gene_list =~ s/${min_id}:${min_dist}:${min_dir}/\*${min_id}:${min_dist}:${min_dir}/;
        }

        $gene_list = $gene_list eq q{.} ? $attribute : "$gene_list; $attribute";

        # print gff line
        print join ("\t",
                     $chr,
                     scalar @range,
                     $feature,
                     $start,
                     $end,
                     $mean,
                     q{.},
                     q{.},
                     $gene_list,
                     "\n");
    }
}



sub gff_filter_by_coord {
    my ($lower_bound, $upper_bound, $data_ref, $overlap) = @_;

    my @filtered;
    for (my $i = 0; $i < @{$data_ref} - 1; $i++) {

        my $start_coord  = $data_ref->[$i][0];
        my $end_coord    = $data_ref->[$i][1];

	if (($end_coord >= $lower_bound && $start_coord <= $upper_bound)
            or ($start_coord <= $lower_bound && $end_coord >= $upper_bound)) {
	    push @filtered, $data_ref->[$i];
	}
	last if ( $start_coord > $upper_bound ); # assumes $data_ref is sorted by coord
    }
    return \@filtered;
}

sub index_gff_annotation {
    my $annotation_file = shift;
    open my $GFFH, '<', $annotation_file or croak "Can't read file: $annotation_file";
    my %annotation = ();
    while (<$GFFH>) {
        next if ($_ =~ m/^#.*$|^\s*$/);
        chomp;
        my %locus = %{gff_read ($_)};

        my ($locus_id) = $locus{attribute} =~ m/ID=([^;]+)/;

        if (!defined $locus_id) {
            ($locus_id, undef) = split /;/, $locus{attribute};
        }

        $locus_id =~ s/["\t\r\n]//g;

        push @{$annotation{$locus{seqname}}},
        [$locus{start}, $locus{end}, $locus_id, $locus{strand}];
    }

    for my $chr (sort {$a cmp $b} keys %annotation) {
        @{$annotation{$chr}}
        = sort {
            $a->[0] <=> $b->[0]
        } @{$annotation{$chr}};
    }

    close $GFFH;
    return \%annotation;
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

 find_loci.pl --- Find overlapping loci in GFF annotation file given GFF probe/window file

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

 find_loci.pl [OPTION]... [FILE]...

 -f,  --ratio-file       gff file with expression, binding sites, etc. data
 -a,  --annotation-file  gff annotation file
 -d,  --distance         distance from center of each probe/window on each side to search
 -p,  --overlap          search whole probe for overlapping genes
 -dp, --distance-overlap search whole probe for overlapping genes and go extra distance beyond bounds
 -w,  --filter           filter loci with length less than the first parameter or more than the second (optional) parameter
 -o,  --output           filename to write results to (default is STDOUT, unless in batch mode)
 -v,  --verbose          output perl's diagnostic and warning messages
 -q,  --quiet            supress perl's diagnostic and warning messages
 -h,  --help             print this information
 -m,  --manual           print the plain old documentation page

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
