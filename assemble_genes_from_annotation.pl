#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

my $gene_id;
my $transcript_id;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'gene-id-field-name|g=s'       => \$gene_id,
    'transcript-id-field-name|t=s' => \$transcript_id,
    'output|o=s' => \$output,
    'verbose|v'  => sub { use diagnostics; },
    'quiet|q'    => sub { no warnings; },
    'help|h'     => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'   => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result and $gene_id and $transcript_id;


if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my %genes = ();

while (<>) {
    chomp;
    next if m/^#.*$|^\s*$/;
    my %locus = %{ gff_read ($_) };
    next unless $locus{feature} =~ m/exon/i;

    my ($gene_id)
    = $locus{attribute} =~ m/.*
                             $gene_id
                             [\s=]+
                             "+
                             ([^;"\t]+)
                             "+
                            /x;

    my ($transcript_id) 
    = $locus{attribute} =~ m/.*
                             $transcript_id
                             [\s=]+
                             "+
                             ([^;"\t]+)
                             "+
                            /x;

    croak "\nCouldn't find the gene id based on the id field name you provided.
Either you provided an invalid field name or my expression matching skills suck...
The line in question is $locus{attribute}.\n\n"
    unless $gene_id;

    croak "\nCouldn't find the transcript id based on the id field name you provided.
Either you provided an invalid field name or my expression matching skills suck...
The line in question is $locus{attribute}.\n\n"
    unless $transcript_id;

    my $chr = $locus{seqname};

    if (!exists $genes{$chr}{$gene_id} or exists $genes{$chr}{$gene_id}{transcript}{$transcript_id}) {
	$genes{$chr}{$gene_id}{transcript}{$transcript_id} = undef;
	$genes{$chr}{$gene_id}{strand}     = $locus{strand }  unless exists $genes{$chr}{$gene_id}{strand};
	$genes{$chr}{$gene_id}{start }     = $locus{start  }  unless exists $genes{$chr}{$gene_id}{start } and $locus{start} > $genes{$chr}{$gene_id}{start};
	$genes{$chr}{$gene_id}{end   }     = $locus{end    }  unless exists $genes{$chr}{$gene_id}{end   } and $locus{start} < $genes{$chr}{$gene_id}{end};
    }
}
    
CHR:
for my $chr (sort keys %genes) {

  GENE:
    for my $gene_id (sort { $genes{$chr}{$a}{start} <=>  $genes{$chr}{$b}{start} } keys %{$genes{$chr}}) {

	my ($transcript_id, @more_than_one_transcript) = keys %{$genes{$chr}{$gene_id}{transcript}};

	if (@more_than_one_transcript > 1) {
	    carp "Found more than one transcript ID for gene $gene_id, skipping:\n", join "\n", @more_than_one_transcript;
	    next GENE;
	}

        print join ("\t",
                    $chr,
                    'dz',
                    'gene',
                    $genes{$chr}{$gene_id}{start},
                    $genes{$chr}{$gene_id}{end},
                    q{.},
                    $genes{$chr}{$gene_id}{strand},
                    q{.},
                    "ID=$gene_id;transcript_id=$transcript_id",
                    "\n"
                );
    }
}


sub gff_read {
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, shift);

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

 assemble_genes_from_annotation.pl - Build gene models from gff exons annotations

=head1 SYNOPSIS

 assemble_genes_from_annotation.pl -g gene_id -t transcript_id -o gene_annotation.gff no_gene_annotation.gff

=head1 DESCRIPTION

=head1 OPTIONS

 assemble_genes_from_annotation.pl [OPTION]... [FILE]...

 -g, --gene-id        name immediately preceeding the gene id in the input gff file
 -t, --transcript-id  name immediately preceeding the transcript id in the input gff file
 -o, --output         filename to write results to (defaults to STDOUT)
 -v, --verbose        output perl's diagnostic and warning messages
 -q, --quiet          supress perl's diagnostic and warning messages
 -h, --help           print this information
 -m, --manual         print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: 249 $:
 $Author: psilva $:
 $Date: 2010-01-11 21:24:34 -0800 (Mon, 11 Jan 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/assemble_genes_from_annotation.pl $:
 $Id: assemble_genes_from_annotation.pl 249 2010-01-12 05:24:34Z psilva $:

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
