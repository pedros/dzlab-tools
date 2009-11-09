#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

my $DATA_HANDLE = 'ARGV';
my $locus_id = 'ID';
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'locus-id|i=s' => \$locus_id,
    'output|o=s'  => \$output,
    'verbose|v'   => sub { use diagnostics; },
    'quiet|q'     => sub { no warnings; },
    'help|h'      => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'    => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my %loci = ();

FEATURE:
while (<$DATA_HANDLE>) {
    next if m/^\s*#/;
    chomp;

    my %locus = %{gff_read ($_)};
    
    my ($locus_id) = $locus{attribute} =~ m/$locus_id[=\s]?([^;]+)/;
    
    my ($c, $t);
    if (!defined $locus_id) {
        ($locus_id, $c, $t) = split /;/, $locus{attribute};
        $locus_id ||= q{.};
    }
    else {
        $locus_id =~ s/["\t\r\n]//g;
    }
    
    if (!defined $c or !defined $t) {
        next FEATURE;
    }

    $loci{$locus}->[0] += $c;
    $loci{$locus}->[1] += $t;
}


for my $locus (sort keys %loci) {
    if ($c + $t > 0) {
        my $score = sprintf ("%g", ($c/($c+$t)));
        print join ("\t", $locus, $score, $c, $t), "\n";
    }
}


sub gff_read {
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute)
    = split /\t/, shift;

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

 $Rev: $:
 $Author: $:
 $Date: $:
 $HeadURL: $:
 $Id: $:

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
