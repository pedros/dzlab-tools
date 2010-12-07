#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;disable diagnostics;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Pod::Usage;

# Globals, passed as command line options
my $ratio_file           = q{-};
my $gff_annotation_file;
my $reference_file;
my $count_CG_sites       = 0;
my $new_feature;
my $no_add               = 0;
my $extend_annotation    = 0;
my $gene_id_field_name   = 'ID';
my $min_score            = 0;
my $average_scores       = 0;
my $output               = q{-};
my $verbose              = 0;
my $quiet                = 0;

my @argv = @ARGV;

# Grabs and parses command line options
my $result = GetOptions (
    'ratio-file|f=s'          => \$ratio_file,
    'gff-annotation-file|g=s' => \$gff_annotation_file,
    'reference-file|r=s'      => \$reference_file,
    'count-CG-sites|cg|c'     => \$count_CG_sites,
    'new-feature|nf=s'        => \$new_feature,
    'no-add|n'                => \$no_add,
    'extend-annotation|e'     => \$extend_annotation,
    'gene-id-field-name|i=s'  => \$gene_id_field_name,
    'min-score|s=i'           => \$min_score,
    'average-scores|a'        => \$average_scores,
    'output|o:s'              => \$output,
    'verbose|v'               => sub {enable diagnostics;use warnings;},
    'quiet|q'                 => sub {disable diagnostics;no warnings;},
    'help|usage|h'            => sub {pod2usage(-verbose => 1);},
    'manual|man|m'            => sub {pod2usage(-verbose => 2);}
);

# Check required command line parameters
pod2usage(-verbose => 1) unless @argv;

# redirects STDOUT to file if specified by user
open STDOUT, '>', "$output" or croak "Can't redirect STDOUT to file: $output" if $output ne q{-};

my %annotation = %{index_gff_annotation ($gff_annotation_file, $gene_id_field_name)} if $gff_annotation_file;

my %data;
open my $RATIO, '<', $ratio_file or croak "Can't read file: $ratio_file" if $ratio_file ne q{-};
while (<$RATIO>) {
    next if ($_ =~ m/^#.*$|^\s*$|\.;\.$/);
    chomp;
    s/[\r\n]//g;
    my $chr = (split /\t/, $_)[0];
    $chr =~ tr/A-Z/a-z/;
    push @{$data{$chr}}, $_;
}
close $RATIO;

my %reference = %{ index_fasta ($reference_file, $count_CG_sites) }
if $reference_file and $count_CG_sites;

CHROMOSOME:
for my $chr (sort {$a cmp $b} keys %annotation) {

    my $last_record = 0;

  ANNOTATION:
    for my $start (sort {$a <=> $b} keys %{$annotation{$chr}}) {
 
        my @range = @{ gff_filter_by_coord ($start, $annotation{$chr}{$start}[1], $last_record, \@{$data{$chr}}, $min_score) };

	$last_record = shift @range;
        
        if ($extend_annotation) {

            if (@range) {
                my ($lowest_coord, $highest_coord)
                = ((split /\t/, $range[0])[3], (split /\t/, $range[$#range])[4]);
                
                if ($lowest_coord < $annotation{$chr}{$start}[0]
                    or $highest_coord > $annotation{$chr}{$start}[1]) {

                    $annotation{$chr}{$start}[6]
                    .= "; Target=$annotation{$chr}{$start}[2] $annotation{$chr}{$start}[0] $annotation{$chr}{$start}[1]";

                    my $original_size
                    = $annotation{$chr}{$start}[1] - $annotation{$chr}{$start}[0];

                    $annotation{$chr}{$start}[0]
                    = $lowest_coord if $lowest_coord < $annotation{$chr}{$start}[0];

                    $annotation{$chr}{$start}[1]
                    = $highest_coord if $highest_coord > $annotation{$chr}{$start}[1];

                    $annotation{$chr}{$start}[6]
                    .= '; extension='
                    . (($annotation{$chr}{$start}[1] - $annotation{$chr}{$start}[0]) - $original_size);

                }

            }

            $annotation{$chr}{$start}[5] = 'ext_' . $annotation{$chr}{$start}[5];

            print join ("\t",
                        $chr,
                        $annotation{$chr}{$start}[4],
                        $new_feature||=$annotation{$chr}{$start}[5],
                        $annotation{$chr}{$start}[0],
                        $annotation{$chr}{$start}[1],
                        q{.},
                        $annotation{$chr}{$start}[3],
                        q{.},
                        $annotation{$chr}{$start}[6],
                    ), "\n";

            next ANNOTATION;
        }


        if ($no_add) {

          WINDOW:
            for my $window (@range) {
                my @fields = split /\t/, $window;

                my $attr = q{.};
                $attr = "$gene_id_field_name=$annotation{$chr}{$start}[2]"
                if $annotation{$chr}{$start}[2] ne q{.};

                print join ("\t",
                            $fields[0],
                            'filtered',
                            $new_feature||$fields[2],
                            $fields[3],
                            $fields[4],
                            $fields[5]||q{.},
                            $fields[6]||q{.},
                            q{.},
                            $attr,
                        ), "\n";
            }
        }
        else {

            my ($actual_score, $actual_attribute, $actual_ct_count, $actual_context);

            unless ($average_scores) {

                my ($context, $a_c_count, $a_t_count, $score, $b_c_count, $b_t_count, $a_score, $b_score)
                = add_gff_attribute_range (\@range);

                my $attribute
                = "$gene_id_field_name=$annotation{$chr}{$start}[2];c=$a_c_count;t=$a_t_count";

                my $locus_len
                = $annotation{$chr}{$start}[1] - $annotation{$chr}{$start}[0] + 1;
            
                my $cent_dist
                = int (($annotation{$chr}{$start}->[1] - $annotation{$chr}{$start}->[0]) / 2) + $annotation{$chr}{$start}->[0];
            
                $cent_dist
                = abs ((length $reference{"$chr-seq"}) - $cent_dist) if $reference_file;

                if (@range == 0) {
                    $attribute = "$gene_id_field_name=$annotation{$chr}{$start}[2]";
                    $score     = q{.};
                }
                elsif (defined $b_c_count) {
                    $score  = sprintf ("%g", $score);
                    $a_score = sprintf ("%g", $a_score);
                    $b_score = sprintf ("%g", $b_score);
                    my $act = sprintf ("%g", $a_c_count + $a_t_count);
                    my $bct = sprintf ("%g", $b_c_count + $b_t_count);
                    $attribute = "$gene_id_field_name=$annotation{$chr}{$start}[2];locus_len=$locus_len;cent_dist=$cent_dist;a_score=$a_score;b_score=$b_score;act_score=$act;bct_score=$bct";
                }
                else {$score = sprintf ("%g", $score)}

                $actual_score = $score;
                $actual_attribute = $attribute;
                $actual_context = $context;
                $actual_ct_count = $a_c_count + $a_t_count;
            }
            else {
                if (@range) {
                    $actual_score = average_gff_score_range (\@range);
                }
                else {
                    $actual_score = q{.};
                }
                $actual_attribute = "$gene_id_field_name=$annotation{$chr}{$start}[2]" || q{.};
            }

            my $total_CG_sites;
            if($count_CG_sites) {
                $total_CG_sites = count_sites (
                    substr(
                        $reference{"$chr-seq"},
                        $annotation{$chr}{$start}->[0],
                        $annotation{$chr}{$start}->[1] - $annotation{$chr}{$start}->[0]
                    ),
                    'CG'
                );
                $actual_attribute
                = "$gene_id_field_name=$annotation{$chr}{$start}->[2];total_CG_sites=$total_CG_sites";
                $actual_attribute .= ";total_ct=$actual_ct_count" if $actual_ct_count;
            }

            print join ("\t",
                        $chr,
                        'dz_win',
                        $new_feature||$actual_context||'window',
                        $annotation{$chr}{$start}->[0],
                        $annotation{$chr}{$start}->[1],
                        $actual_score,
                        $annotation{$chr}{$start}->[3],
                        q{.},
                        $actual_attribute,
                    ), "\n";

         }
    }
}


sub count_sites {
    my ($sequence, $type) = @_;

    my @count = $sequence =~ m/$type/g;

    return @count;
}


sub average_gff_score_range {
    my $range = shift;

    my ($score_total, $score_count) = (0, 0);

    foreach my $k (@{$range}) {
        my $current_rec = gff_read ($k);

        next if $current_rec->{score} !~ m/\d/;

        $score_total += $current_rec->{score};
        $score_count++;
    }

    return ($score_count ? $score_total / $score_count : 'NaN');
}

sub add_gff_attribute_range {
    my $range = shift;
    my ($a_c_count, $a_t_count, $b_c_count, $b_t_count, $score, $a_score, $b_score, $context)
    = (0, 0, 0, 0, 0, 0, 0, q{.});

    my ($score_count, $score_total) = (0, 0);

    my @attr_fields = ();
    foreach my $k (@{$range}) {
        my %current_rec = %{gff_read ($k)};

        next if $current_rec{attribute} eq q{.};

        @attr_fields = split(/;/, $current_rec{'attribute'});
        my ($p_value, $a_c_tmp, $a_t_tmp, $b_c_tmp, $b_t_tmp);

        if (@attr_fields == 2) {
            ($a_c_tmp, $a_t_tmp) = @attr_fields;
        }
        elsif (@attr_fields == 5) {
            ($p_value, $a_c_tmp, $a_t_tmp, $b_c_tmp, $b_t_tmp) = @attr_fields;
        }
        else {
            croak "Wrong number of attribute fields: must be 2 (for single c gff files) or 5 (for inter-tissue comparison files.";
        }

        ($a_c_tmp) = $a_c_tmp =~ m/(\d+)/;
        ($a_t_tmp) = $a_t_tmp =~ m/(\d+)/;
        $a_c_count += $a_c_tmp;
        $a_t_count += $a_t_tmp;

        if (@attr_fields == 5) {
            ($b_c_tmp) = $b_c_tmp =~ m/(\d+)/;
            ($b_t_tmp) = $b_t_tmp =~ m/(\d+)/;
            $b_c_count += $b_c_tmp;
            $b_t_count += $b_t_tmp;
        }

        croak (
            "Found different site context when merging ranges:
             last site was $context and new site is $current_rec{feature}."
        ) if $context ne q{.} and $context ne $current_rec{feature};

        $context = $current_rec{feature};
    }

    my @scores = ();
    if (@attr_fields < 5 and !$score_total) {
        $score = $a_c_count / ($a_c_count + $a_t_count) if $a_c_count + $a_t_count != 0;
        @scores = ($a_c_count, $a_t_count, $score);
    }
    else {
        $a_score = $a_c_count / ($a_c_count + $a_t_count) if $a_c_count + $a_t_count != 0;
        $b_score = $b_c_count / ($b_c_count + $b_t_count) if @attr_fields == 5 and $b_c_count + $b_t_count != 0;
        $score = $a_score - $b_score;
        @scores = ($a_c_count, $a_t_count, $score, $b_c_count, $b_t_count, $a_score, $b_score);
    }
    return ($context, @scores);
}


sub gff_filter_by_coord {
    my ($lower_bound, $upper_bound, $last_index_seen, $data_ref, $min_score) = @_;

    my @filtered;
  WINDOW:
    for (my $i = $last_index_seen; $i < @{$data_ref}; $i++) {

        my ($start_coord, $end_coord, $score)
        = (split /\t/, $data_ref->[$i])[3, 4, 5];

	if ($end_coord >= $lower_bound && $start_coord <= $upper_bound) {
          
            next WINDOW if $min_score and ($score !~ m/\d/ or $score < $min_score);
  
	    push @filtered, $data_ref->[$i];
	    $last_index_seen = $i;
	}

	last if ( $start_coord > $upper_bound );

    }
    unshift @filtered, $last_index_seen;
    return \@filtered;
}


sub index_gff_annotation {
    my ($annotation_file, $gene_id_field_name) = @_;

    open my $GFFH, '<', $annotation_file or croak "Can't read file: $annotation_file";
    my %annotation = ();
    while (<$GFFH>) {
        next if ($_ =~ m/^#.*$|^\s*$/);
        chomp;
        s/[\r\n]//g;
        my %locus = %{gff_read ($_)};

        my ($locus_id) = $locus{attribute} =~ m/$gene_id_field_name[=\s]?([^;]+)/;

        if (!defined $locus_id) {
            ($locus_id, undef) = split /;/, $locus{attribute};
            $locus_id ||= q{.};
        }
        else {
            $locus_id =~ s/["\t\r\n]//g;
        }

        $annotation{$locus{seqname}}{$locus{start}}
        = [$locus{start}, $locus{end}, $locus_id, $locus{strand}, $locus{source}, $locus{feature}, $locus{attribute}];
    }
    close $GFFH;
    return \%annotation;
}


sub gff_read {
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute)
    = split(/\t/, shift);

    my %rec = (
	'seqname'   => lc $seqname,
	'source'    => $source,
	'feature'   => $feature,
	'start'     => $start,
	'end'       => $end,
	'score'     => $score,
	'strand'    => $strand,
	'frame'     => $frame,
	'attribute' => $attribute
	);
    return \%rec;
}

sub index_fasta {
    my ($reference_file, $count_CG_Sites) = @_;

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

        my ($separator) = $line =~ m/(N{1000,})/i; # looks for regions with more than 1000 Ns
        if ($separator) {
            my $sep_start = (index $line, $separator) + 1;
            my $sep_end   = length $separator;
            my $sep_cen   = ($sep_end / 2) + $sep_start;
            $reference{$dsc[$j]} = $sep_cen;
#            $line =~ tr/ACGTacgt/TGCAtgca/ if $count_CG_sites;
#            $reference{"$dsc[$j]-rc"} = reverse $line if $count_CG_sites;
        }
        else {
            print STDERR "No centrometer region found for $dsc[$j].\n";
        }
        $reference{"$dsc[$j]-seq"} = $line if $count_CG_sites;
    }
    return \%reference;
}


__END__

=head1 NAME

 window_by_annotation.pl - Window input GFF file to GFF annotated loci

=head1 VERSION

 $Rev: 249 $:
 $Author: psilva $:
 $Date: 2010-01-11 21:24:34 -0800 (Mon, 11 Jan 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/window_by_annotation.pl $:
 $Id: window_by_annotation.pl 249 2010-01-12 05:24:34Z psilva $:

=head1 USAGE

 # window input.gff to loci in annotation.gff averaging methylation information present in attribute
 # and reporting number of CG sites per locus present in the reference genome
 window_by_annotation.pl -f input.gff -g annotation.gff -r reference.fasta -c

 # report which locus each input window falls into. extract loci IDs using 'transcript_id' GFF attribute tag
 window_by_annotation.pl -f input.gff -g annotation.gff -n -i transcript_id

 # extend annotation.gff using the windows in input.gff that overlap each locus boundaries
 # extract loci IDs using 'ID' GFF attribute tag
 window_by_annotation.pl -f input.gff -g annotation.gff -e -i ID

=head1 OPTIONS

 window_by_annotation.pl [OPTION]... [FILE]...

 -f,  --ratio-file          GFF data input file
 -g.  --gff-annotation-file GFF annotation file
 -r,  --reference-file      FASTA reference file
 -c,  --count-CG-sites      compute number of CG sites per locus
 -nf, --new-feature         substitute feature ID
 -n,  --no-add              don't try to sum up numeric attributes
 -a,  --average-scores      disregard attributes fields; only average scores
 -e,  --extend-annotation   'reverse' windowing: extend annotations using overlapping GFF data
 -i,  --gene-id-field-name  GFF attribute locus interest tag [ID]
 -o, --output      filename to write results to (default is STDOUT, unless in batch mode)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 DESCRIPTION

 Takes a GFF-formatted input file with methylation information on the attribute field
 in the form c=? and t=? *or* a generic GFF data file.
 Given an annotation file also in GFF format, windows all records in input per each locus.
 It can also extend the annotation file using the overlapping GFF data information

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 LICENSE AND COPYRIGHT

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
