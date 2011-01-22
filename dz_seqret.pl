#!/usr/bin/env perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell
# $Id: dz_seqret.pl 249 2010-01-12 05:24:34Z psilva $
# -*-Perl-*- mode (for emacs)

=head1 NAME

bp_seqret - bioperl implementation of sequence fetch from local db (like EMBOSS seqret)

=head1 USAGE

bp_seqret [-f/--format outputformat] [-o/--out/--outfile outfile] [-d/--db dbname] [-i/--id/-s/--seqname seqname1]

Example usage:

   bp_seqret -f fasta -db db.fa -i seq1 -i seq2 > output.fas
   bp_seqret db.fa:seq1 output.fas
   bp_seqret db.fa:seq1 -o output.fas
   bp_seqret -db db.fa -o output.fas seq1 seq2 seq3
   bp_seqret -db db.fa seq1 seq2 seq3 output.fas
   bp_seqret -db db.fa seq1 seq2 seq3 - > output.fas  

The DB is expected to be a Fasta formatted sequence file with multiple
sequences.

Output format is Fasta by default.

If no output filename is provided then output is written to STDOUT.
Providing '-' as the output filename will accomplish the same thing.


=head1 AUTHOR

Jason Stajich jason_AT_bioperl-dot-org

=cut

use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Getopt::Long;

my $dbname;
my $format = 'fasta';
my $outfile;
my ($start,$end);
my $gff;
my $feature;
GetOptions(
	   'f|format:s'   => \$format,
	   'o|out|outfile:s' => \$outfile,
	   's|sbegin|begin|start:s'  => \$start,
	   'e|send|end|stop:s'       => \$end,
	   'd|db|dbname:s'   => \$dbname,
	   'i|id|seqname:s'  => \$gff,
           'feature=s'       => \$feature
);

if( ! $dbname ) {
    die "need a dbname\n" unless @ARGV;
    $dbname = shift @ARGV;	
    # if( $dbname =~ s/^([^:]+):// ) {
    #     push @names, $dbname;
    #     $dbname = $1;
    # }				
}

my $db = Bio::DB::Fasta->new($dbname, -glob => "*.{fa,fas,fsa,fasta,pep,aa,seq,cds,peps}");
if( ! $outfile ) {
    $outfile = pop @ARGV;
}
my $out;
if( $outfile ) {
    $out = Bio::SeqIO->new(-format => $format,
			   -file   => ">$outfile");
} else {
    $out = Bio::SeqIO->new(-format => $format);
}

open my $GFF, '<', $gff, or die "Can't open $gff: $!";
while (<$GFF>) {
    next if m/^\s*#|^\s*$/;
    chomp;

    my ($chr, $feat, $start, $end, $strand)
    = (split /\t/, $_)[0, 2, 3, 4, 6];

    next if $feature && $feat ne $feature;

    if ($strand eq q{-}) {
        $start ^= $end; $end ^= $start; $start ^= $end;
    }

    my $seq;
    if( $start || $end ) {
	$seq = $db->seq($chr, $start => $end);
    } else { 
	$seq = $db->seq($chr);
    }
    if( $seq ) { 
	my ($id,$desc) = split(/\s+/,$db->header($chr),2);
	if( $start && $end ) { 
	    $id = sprintf("%s_%d-%d",$id,$start || 0,$end || 0);
	}
	
	$out->write_seq(Bio::PrimarySeq->new(-display_id => $id,
					     -description => $desc,
					     -seq => $seq));
    } else {
	warn("$chr not found\n");
    }
}

