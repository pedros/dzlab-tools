#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;disable diagnostics;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Pod::Usage;

# Globals, passed as command line options
my $ratio_file = q{-};
my $gff_annotation_file;
my $annotation_file;
my $reference_file;
my $no_add  = 0;
my $output  = q{-};
my $verbose = 0;
my $quiet   = 0;

my @argv = @ARGV;

# Grabs and parses command line options
my $result = GetOptions (
    'ratio-file|f=s'      => \$ratio_file,
    'gff-annotation-file|g=s' => \$gff_annotation_file,
    'annotation-file|a=s' => \$annotation_file,
    'reference-file|r=s'  => \$reference_file,
    'no-add|n'            => \$no_add,
    'output|o:s'          => \$output,
    'verbose|v'           => sub {enable diagnostics;use warnings;},
    'quiet|q'             => sub {disable diagnostics;no warnings;},
    'help|usage|h'        => sub {pod2usage(-verbose => 1);},
    'manual|man|m'        => sub {pod2usage(-verbose => 2);}
);

# Check required command line parameters
pod2usage(-verbose => 1) unless @argv;

# redirects STDOUT to file if specified by user
open STDOUT, '>', "$output" or croak "Can't redirect STDOUT to file: $output" if $output ne q{-};

my %annotation = ();
%annotation = %{index_gff_annotation ($gff_annotation_file)} if $gff_annotation_file;
%annotation = %{index_generic_annotation ($annotation_file)} if $annotation_file;

my %data;
open my $RATIO, '<', $ratio_file or croak "Can't read file: $ratio_file" if $ratio_file ne q{-};
while (<$RATIO>) {
    chomp;
    next if ($_ =~ m/^#.*$|^\s*$|\.;\.$/);
    my $chr = (split /\t/, $_)[0];
    $chr =~ tr/A-Z/a-z/;
    push @{$data{$chr}}, $_;
}
close $RATIO;

# prints out header fields that contain gff v3 header, generating program, time, and field names
# gff_print_header ($0, @argv);

my %reference = %{ index_fasta ($reference_file) };

for my $chr (sort {$a cmp $b} keys %annotation) {

    my $last_record = 0;

    for my $start (sort {$a <=> $b} keys %{$annotation{$chr}}) {

        my @range = @{ gff_filter_by_coord ($start, $annotation{$chr}{$start}[1], $last_record, \@{$data{$chr}}) };

	$last_record = shift @range;

        if ($no_add) {
            for my $window (@range) {
                my @fields = split /\t/, $window;
                my ($p, $ac, $at, $bc, $bt) = (split /;/, $fields[-1]);
                ($p)  = $p  =~ m/(\d+)/;
                ($ac) = $ac =~ m/(\d+)/;
                ($at) = $at =~ m/(\d+)/;
                ($bc) = $bc =~ m/(\d+)/;
                ($bt) = $bt =~ m/(\d+)/;

                my ($as, $bs) = (0, 0);
                $as = $ac / ($ac + $at) if $ac + $at != 0;
                $bs = $bc / ($bc + $bt) if $bc + $bt != 0;

                my $sc   = sprintf ("%6f", $as - $bs);
                $as = sprintf ("%g", $as);
                $bs = sprintf ("%g", $bs);
                my $a_ct  = sprintf ("%g", $ac + $at);
                my $b_ct  = sprintf ("%g", $bc + $bt);

                my $attr;
                if ($reference_file) {
                    my $locus_len = $annotation{$chr}{$start}[1] - $annotation{$chr}{$start}[0] + 1;
                    my $cent_dist = abs ($reference{$fields[0]} - ( int (($fields[4] - $fields[3]) / 2) + $fields[3]));
                    $attr = "ID=$annotation{$chr}{$start}[2]\t$locus_len\t$cent_dist\t$as\t$bs\t$a_ct\t$b_ct";
                }
                else {
                    $attr = "$as\t$bs\t$a_ct\t$b_ct";
                }

                print join ("\t",
                            $fields[0],
                            'filtered',
                            $fields[2],
                            $fields[3],
                            $fields[4],
                            $sc,
                            q{.},
                            q{.},
                            $attr,
                        ), "\n";
            }
        }
        else {
            my ($context, $a_c_count, $a_t_count, $score, $b_c_count, $b_t_count, $a_score, $b_score)
            = add_gff_attribute_range (\@range);

            my $attribute = "c=$a_c_count;t=$a_t_count";

            my $locus_len = $annotation{$chr}{$start}[1] - $annotation{$chr}{$start}[0] + 1;
            my $cent_dist = int (($annotation{$chr}{$start}->[1] - $annotation{$chr}{$start}->[0]) / 2) + $annotation{$chr}{$start}->[0];
            $cent_dist = abs ($reference{$chr} - $cent_dist) if $reference_file;

            if (@range == 0) {
                $attribute = q{.};
                $score     = q{.};
            }
            elsif (defined $b_c_count) {
                $score  = sprintf ("%g", $score);
                $a_score = sprintf ("%g", $a_score);
                $b_score = sprintf ("%g", $b_score);
                my $act = sprintf ("%g", $a_c_count + $a_t_count);
                my $bct = sprintf ("%g", $b_c_count + $b_t_count);
                $attribute = "ID=$annotation{$chr}{$start}[2]\t$locus_len\t$cent_dist\t$a_score\t$b_score\t$act\t$bct";
            }
            else {$score = sprintf ("%g", $score)}

            print join ("\t",
                        $chr,
                        'window',
                        $context,
                        $annotation{$chr}{$start}->[0],
                        $annotation{$chr}{$start}->[1],
                        $score,
                        q{.},
                        q{.},
                        $a_c_count + $a_t_count,#$attribute,
                        "ID=$annotation{$chr}{$start}->[2]",
                    ), "\n";

#             print join ("\t",
#                         $annotation{$chr}{$start}->[2],
#                         $annotation{$chr}{$start}->[3],
#                         $annotation{$chr}{$start}->[0],
#                         $annotation{$chr}{$start}->[1],
#                         $a_c_count,
#                         $a_t_count,
#                         $score,
#                     ), "\n";
         }
    }
}


sub add_gff_attribute_range {
    my $range = shift;
    my ($a_c_count, $a_t_count, $b_c_count, $b_t_count, $score, $a_score, $b_score, $context) = (0, 0, 0, 0, 0, 0, 0, q{.});

    my @attr_fields = ();
    foreach my $k (@{$range}) {
        my %current_rec = %{gff_read ($k)};

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
    if (@attr_fields == 2) {
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
    my ($lower_bound, $upper_bound, $last_index_seen, $data_ref) = @_;

    my @filtered;
    for (my $i = $last_index_seen; $i < @{$data_ref}; $i++) {

        my $start_coord = (split /\t/, $data_ref->[$i])[3];

	if ($start_coord >= $lower_bound && $start_coord <= $upper_bound) {
	    push @filtered, $data_ref->[$i];
	    $last_index_seen = $i;
	}

	last if ( $start_coord > $upper_bound );

    }
    unshift @filtered, $last_index_seen;
    return \@filtered;
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

sub annon_read {
    my ($locus_id, $comment, $start, $end, $length) = split(/\t/, shift);

    $locus_id =~ m/^at(\d?)/i;

    my $seqname;
    $seqname = 'chr' . $1 if $1;
    $seqname = $comment unless $1;


    my %rec = (
        'seqname'  => $seqname,
        'locus_id' => $locus_id,
        'comment'  => $comment,
	'start'    => $start,
	'end'      => $end,
        'length'   => $length,
    );

    return \%rec;
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

        $annotation{$locus{seqname}}{$locus{start}} = [$locus{start}, $locus{end}, $locus_id];
    }
    close $GFFH;
    return \%annotation;
}


sub index_generic_annotation {
    my $annonfile = shift;
    open my $ANNON, '<', $annonfile or croak "Can't open $annonfile.";

    my %annotation = ();

    while (<$ANNON>) {
        next if ($_ =~ m/^#.*$|^\s*$/);
        chomp;
        my %locus = %{annon_read ($_)};

        $annotation{$locus{seqname}}{$locus{start}} = [$locus{start}, $locus{end}, $locus{locus_id}, $locus{comment}];
    }
    close $ANNON;
    return \%annotation;
}


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

        my ($separator) = $line =~ m/(N{1000,})/i; # looks for regions with more than 1000 Ns
        if ($separator) {
            my $sep_start = (index $line, $separator) + 1;
            my $sep_end   = length $separator;
            my $sep_cen   = ($sep_end / 2) + $sep_start;
            $reference{$dsc[$j]} = $sep_cen;
        }
        else {
            print STDERR "No centrometer region found for $dsc[$j].\n";
        }
    }
    return \%reference;
}

sub gff_print_header {
    my @call_args = @_;
    print "##gff-version 3\n";
    print join(' ',
	       '#',
	       @call_args,
	       "\n");
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime (time);
    printf "# %4d-%02d-%02d %02d:%02d:%02d\n", $year+1900, $mon+1, $mday, $hour, $min, $sec;

    print join("\t",
	       '# SEQNAME',
	       'SOURCE',
	       'FEATURE',
	       'START',
	       'END',
	       'SCORE',
	       'STRAND',
	       'FRAME',
               'LOCUSID',
               'LOCUSLEN',
               'WINCENTDIST',
	       'AC/T',
               'BC/T',
               'AC+T',
               'BC+T',
	), "\n";
    return 0;
}

