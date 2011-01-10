#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use version; our $VERSION = qv('0.0.1');

use FindBin;
use lib "$FindBin::Bin/DZLab-Tools/lib";
use DZLab::Tools::RunUtils;

use Pod::Usage;
use Getopt::Long;

my $help;
my $verbose;
my $annotation;
my $result = GetOptions (
    "verbose" => \$verbose,
    "help"    => \$help,
    "annotation|a=s" => \$annotation,
);
pod2usage(-verbose => 1) if (!$annotation || $help);  

my %genes;
open my $ANNOTATION, q{<}, $ARGV{annotation} or die "Can't open $ARGV{annotation}: $!";

LOAD_CANONICAL_GENE_LIST:
while (<$ANNOTATION>) {
    chomp;
    my ($id) = $_ =~ m/ID=([^;]+)/;
    $genes{$id} = 0;
}
close $ANNOTATION or die "Can't close $ARGV{annotation}: $!";

my @genes = sort keys %genes;

SCORES_FILE:
for my $scores (@ARGV) {
    my %file_genes;
    open my $SCORES, q{<}, $scores or die "Can't open $scores: $!";

  SCORE:
    while (<$SCORES>) {
        chomp;
        my ($id, $sc, $alt) = split /\t/; $alt ||= '';

        $file_genes{$id} = join "\t", $sc, $alt;
    }
    close $SCORES or die "Can't close $scores: $!";

    open my $OUTFILE, q{>}, "$scores.allgenes" or die "Can't open $scores: $!";
  CANONICAL_GENE:
    for (@genes) {
        if (exists $file_genes{$_}) {
            print $OUTFILE "$_\t$file_genes{$_}\n";
        }
        else {
            print $OUTFILE "$_\t$genes{$_}\n";
        }
    }
    close $OUTFILE; 
}

__DATA__


__END__

=head1 NAME

 ensure_all_genes_present.pl - Ensure all genes present

=head1 SYNOPSIS

 ensure_all_genes_present.pl -a annotation.gff input_file1.gff inputfile2.gff ...

=head1 DESCRIPTION

For each input file, output entire content PLUS a line for every gene which 
is present in the annotation file but not in the input file. Output file is 
named input_fileN.gff.allgenes

=head1 OPTIONS

     --annotation      GFF Annotation file to compare the input file genes against.  
      -a               Each line in annotation file is parsed for and ID=...; field.

     --verbose         print increasingly verbose error messages
     --quiet           print no diagnostic or warning messages
     --version         print current version
     --license         print author's contact and copyright information
     --help            print this information
     --manual          print the plain old documentation page

=head1 VERSION

 0.0.2

=head1 REVISION

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
