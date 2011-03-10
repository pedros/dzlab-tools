#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser;
use Log::Log4perl qw/:easy/;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) if !$opt_input;

Log::Log4perl->easy_init({ 
        level    => $INFO,
        layout   => "%d{HH:mm:ss} %p> (%L) %M - %m%n", 
        #file     => ">log4perl.log",
    });

my $gene_tag = $opt_gene_tag;
my $exon_tag = $opt_exon_tag;
my $single   = $opt_single;
my $in_file  = $opt_input;
my $out_file = $opt_output;

my %genes = (); 
# {
#   seq  => { 
#     locus => { 
#           gff => $gff, 
#           isoform => {
#             isoform1 => [$gff, ...], 
#             isoform2 => [$gff, ...] 
#           }
#     }
#   }
# }

open my $out, '>', $out_file;

my $parser = GFF::Parser->new(file => $in_file);

my $counter=0;
while (my $gff = $parser->next){
    if (++$counter % 5000 == 0) { INFO("$counter"); }

    if ($gff->feature eq 'gene'){
        if (my $gene_locus = $gff->get_attribute($gene_tag)){
            $genes{$gff->sequence}{$gene_locus}{gff} = $gff;
        }
    } 
    elsif ($gff->feature eq 'exon'){
        if (my $exon_locus = $gff->get_attribute($exon_tag)){
            my ($gene_locus, $id) = $gff->parse_locus($exon_tag);
            if ($gene_locus && (!$single || $id == 1)){
                push @{$genes{$gff->sequence}{$gene_locus}{isoform}{$exon_locus}}, $gff;
            }
        }
    } 
}

for my $seq (sort {$a cmp $b} keys %genes) {
    for my $locus (sort {
            $genes{$seq}{$a}{gff}->start <=> $genes{$seq}{$b}{gff}->start
        } grep {
            exists $genes{$seq}{$_}{gff} && exists $genes{$seq}{$_}{isoform}
        } keys %{$genes{$seq}}) {

        # {gff => $gff, isoform => {isoform1 => $gff, ..}}
        my $locus_hash = $genes{$seq}{$locus}; 
        my $gene_gff   = $locus_hash->{gff};

        my $start      = $gene_gff->start;
        my $end        = $gene_gff->end;

        DEBUG($gene_gff->to_string);

        for my $exon_locus (sort keys %{$locus_hash->{isoform}}){
            my @exons = sort {$a->start <=> $b->start} @{$locus_hash->{isoform}{$exon_locus}};
            for my $exon (@exons){
                DEBUG($exon->to_string);
            }
            for my $exon (@exons){
                if ($start < $exon->start){
                    #DEBUG( 
                    say $out (
                        join "\t", 
                        $gene_gff->sequence, 
                        $gene_gff->source, 
                        'intron', 
                        $start, 
                        ($exon->start-1),
                        q{.},
                        $gene_gff->frame // q{.},
                        q{.},
                        $exon->attribute_string,
                    );
                }
                $start = $exon->end + 1;
            }
        }
    }
}

=head1 NAME

make_intron_gff.pl - Create intron file from exon and gene annotation file.

=head1 SYNOPSIS

 make_intron_gff.pl --input TAIR8_gmod.gff --output introns.gff --single

=head1 OPTIONS

=over

=item --gene-tag <tag> | -g <tag>

Gene's locus key. (default: ID)

=for Euclid 
    tag.default: 'ID'

=item --exon-tag <tag> | -g <tag>

Exon's locus key. (default: Parent)

=for Euclid 
    tag.default: 'Parent'

=item  --single

grab only one isoform per gene.

=item --input <file>

Input GFF File with exons and gffs.

=for Euclid
    file.type:  readable

=item --output <file>

Output file (stdout by default)

=for Euclid
    file.default:  '-'

=back

=cut

