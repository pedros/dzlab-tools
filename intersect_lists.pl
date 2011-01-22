#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

my $output;
my $annotation;
my $lists_tag = 'gene';
my $genes_tag = 'ID';

# Grabs and parses command line options
my $result = GetOptions (
    'annotation|a=s' => \$annotation,
    'lists-tag|lt=s' => \$lists_tag,
    'genes-tag|gt=s' => \$genes_tag,
    'output|o=s'     => \$output,
    'verbose|v'      => sub { use diagnostics; },
    'quiet|q'        => sub { no warnings; },
    'help|h'         => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'       => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $result and @ARGV;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my $genes_regex = qr/$genes_tag[=\s]?([^;\s]+).*/;
my $lists_regex = qr/$lists_tag[=\s]?.*\*([^;:\s]+).*/;

my $genes         = index_genes ($annotation, $genes_regex);
my $lists         = {};
my $sets          = {};
my $intersections = {};

while (my $list = shift @ARGV) {

    $lists->{$list} = index_genes ($list, $lists_regex);

    for my $gene (sort keys %{$genes}) {
        if (exists $lists->{$list}->{$gene}) {
            push @{$genes->{$gene}{list}}, $list;
            $sets->{$list}++;
        }
    }
}

my $output_files = {};

for my $gene (sort keys %{$genes}) {
    next unless @{$genes->{$gene}{list}};
    my $combo = join q{,}, sort @{$genes->{$gene}{list}};
    $intersections->{$combo}++;

    $output_files->{$combo}
    || open $output_files->{$combo}, '>', "$combo.gene-list.gff"
    || croak "Can't write to $combo: $!";
    
    print {$output_files->{$combo}} "$genes->{$gene}{gff}\n";
}

print q{=} x 40, "\n";
print "Genes per list:\n";
print q{=} x 40, "\n";

while (my ($list, $count) = each %$sets) {
    print $list, "\t", $count, "\n";
}

print "\n";
print q{=} x 40, "\n";
print "Genes per intersection:\n";
print q{=} x 40, "\n";

while (my ($list, $count) = each %$intersections) {
    print $list, "\t", $count, "\n";
}




sub index_genes {
    my ($genes_file, $regex) = @_;

    my $genes = {};

    open my $GENES, '<', $genes_file or croak "Can't open $genes_file: $!";
    while (<$GENES>) {
        next if m/^\s*#/;
        s/[\r\n]//g;
        my @fields = split /\t/;

        $fields[-1] =~ s/$regex/$1/;

        $genes->{$fields[-1]}{list} = [];
        $genes->{$fields[-1]}{gff} = $_;
    }

    return $genes;
}



__END__


=head1 NAME

 intersect_lists.pl - Given a GFF annotation file and a list of GFF data files, determine which data files contain which loci

=head1 SYNOPSIS

 intersect_lists.pl --genes TAIR8_genes.gff list_1.gff list_2.gff ... list_n.gff

 intersect_lists.pl -g dmel_genes.gff list*gff --output dmel_lists_intersection.dat

=head1 DESCRIPTION

 Given a GFF annotation file and a list of GFF data files, determines the intersection of each data file and the loci in the annotation.

 It requires a valid GFF version 3 file format, in which the last (9th) field--the attribute--may contain multiple TAG=VALUE pairs delimited by a semi-collon ':'.

 It assumes by default that, for the GFF annotation file, that the field to compare is TAGged by 'ID=', but this can be changed via the --genes-tag option
 (eg. --genes-tag 'transcript_id').

 It assumes by default that the GFF data files have attributes delimited by 'gene=', that there may be multiple values per TAG delimited by a comma ',',
 and that in that case the appropriate value will be determined by its having a '*' between the tag and the value (eg. gene=VALUE1,*VALUE2 will return VALUE2).

=head1 OPTIONS

 intersect_lists.pl [OPTION]... [FILES]...

 -a,  --annotation  filename of loci GFF file
 -lt, --lists-tag   GFF attribute identifier tag (default: ID)
 -gt, --genes-tag   GFF attribute identifier tag (default: gene)
 -o,  --output      filename to write results to (defaults to STDOUT)
 -v,  --verbose     output perl's diagnostic and warning messages
 -q,  --quiet       supress perl's diagnostic and warning messages
 -h,  --help        print this information
 -m,  --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: 381 $:
 $Author: psilva $:
 $Date: 2010-07-08 15:53:27 -0700 (Thu, 08 Jul 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/intersect_lists.pl $:
 $Id: intersect_lists.pl 381 2010-07-08 22:53:27Z psilva $:

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
