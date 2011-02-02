#!/usr/bin/env perl
# ___UNDOCUMENTED___

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my @lists = ();
my $gene_methylation_file;
my $gene_id_field_name = 'ID';
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'annotation-file|a=s' => \$gene_methylation_file,
    'lists|l=s{,}'        => \@lists,
    'output|o=s'          => \$output,
    'verbose|v'           => sub { use diagnostics; },
    'quiet|q'             => sub { no warnings; },
    'help|h'              => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'            => sub { pod2usage ( -verbose => 2 ); }
);

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

my $genes = index_gff_annotation ($gene_methylation_file, $gene_id_field_name);

my ($percentile, $missed_gene_count) = (1, 0);

for my $list (@lists) {

    my ($total_length, $total_score, $total_cgsite, $total_ctsite,
        $total_cg_adjusted_score, $total_genes, $total_c, $total_t)
    = (0, 0, 0, 0, 0, 0, 0, 0);

    open my $LIST, '<', $list or croak "Can't open $list";
  GENE:
    while (my $gene = <$LIST>) {

	chomp $gene;
	my ($gene_id, $freq, $alt);

        my @fields = split /\t/, $gene;

        if (@fields < 9) {
            ($gene_id, $freq, $alt) = @fields;
        }
        else {
            ($gene_id, $freq, $alt) = @fields[8, 5, 0];
            $gene_id =~ s/^.*$gene_id_field_name=([^;]+).*$/$1/;
        }

        unless (exists $genes->{$gene_id}) {
            # print STDERR "Can't find ID $gene_id in ", (split m{/}, $gene_methylation_file)[-1], "\n";
            $missed_gene_count++;
            next GENE;
        }

	if ($genes->{$gene_id}->[1] eq 'NaN') {
	    next GENE;
	}

	$total_length += $genes->{$gene_id}->[0];
	$total_score  += $genes->{$gene_id}->[1];
	$total_cgsite += $genes->{$gene_id}->[2];
	$total_ctsite += $genes->{$gene_id}->[3] if defined $genes->{$gene_id}->[3];
        $total_cg_adjusted_score += $total_score * $total_cgsite;

        $total_genes++;

	if ($total_ctsite) {
	    $total_c += ($genes->{$gene_id}->[1] * $genes->{$gene_id}->[3]);
	    $total_t += $genes->{$gene_id}->[3] - ($genes->{$gene_id}->[1] * $genes->{$gene_id}->[3]);
	}
    }
    close $LIST;

    my $arithmetic_mean  = ($total_genes ? ($total_score / $total_genes) : 'NaN');
    my $fractional_meth  = ($total_ctsite and ($total_c or $total_t)
			    ? ($total_c / ($total_c + $total_t))
			    : 'NaN');
    my $cg_adjusted_mean = ($total_cgsite ? ($total_cg_adjusted_score / $total_cgsite) : 'NaN');

    print join ("\t",
                $percentile++,
                $arithmetic_mean,
                $fractional_meth,
                $cg_adjusted_mean,
                $total_length,
                $total_genes,
            ), "\n";

    if ($missed_gene_count) {
	# print STDERR "Couldn't find $missed_gene_count out of " . ($missed_gene_count + $total_genes) . " IDs in ",
	# (split m{/}, $gene_methylation_file)[-1], "\n";
	$missed_gene_count = 0;
    }
}

## done


sub index_gff_annotation {
    my ($annotation_file, $gene_id_field_name) = @_;

    my %annotation = ();

    open my $GFFH, '<', $annotation_file or croak "Can't read file: $annotation_file";
    while (<$GFFH>) {
        next if ($_ =~ m/^#.*$|^\s*$/);
        chomp;
        my %locus = %{gff_read ($_)};

        next if $locus{feature} eq q{.};

        my ($locus_id) = $locus{attribute} =~ m/.* $gene_id_field_name [\s=] "? (\w*\d*) "?/x;
        my ($CG_sites) = $locus{attribute} =~ m/.* total_CG_sites [=\s] "? ([^;]+) "?/x;
        my ($ct_sites) = $locus{attribute} =~ m/.* total_ct [=\s] "? ([^;]+) "?/x;

        unless (defined $locus_id) {
            ($locus_id, undef, undef) = split /;/, $locus{attribute};
            $locus_id =~ s/["\t\r\n]//g;
        }

        # next unless $CG_sites and $ct_sites;

        $annotation{$locus_id}
        = [ ($locus{end} - $locus{start} + 1), $locus{score}, $CG_sites, $ct_sites ];
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

 $Rev: 249 $:
 $Author: psilva $:
 $Date: 2010-01-11 21:24:34 -0800 (Mon, 11 Jan 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/average-partitions.pl $:
 $Id: average-partitions.pl 249 2010-01-12 05:24:34Z psilva $:

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
