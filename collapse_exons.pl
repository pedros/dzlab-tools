#!/usr/bin/env perl
# ___UNDOCUMENTED___

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

my $DATA_HANDLE = 'ARGV';
my $genes;
my $exons;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'genes|g=s'    => \$genes,
    'exons|e=s'    => \$exons,
    'output|o=s'   => \$output,
    'verbose|v'    => sub { use diagnostics; },
    'quiet|q'      => sub { no warnings; },
    'help|h'       => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'     => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $result;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my %exons = %{ index_gff ($exons) };
my %genes = %{ index_gff ($genes) };


for my $locus (sort keys %exons) {

    next unless exists $genes{$locus};

    my ($ec, $et) = @{$exons{$locus}};
    my ($gc, $gt) = @{$genes{$locus}};
    my ($escore, $gscore);

    # if ($ec + $et > 0) {
    #     $escore = sprintf ("%g", ($ec/($ec+$et)));
    #     print join ("\t", $locus, $escore, $ec, $et), "\n";
    # }

    if ($gc + $gt > 0) {
        $gscore = sprintf ("%g", ($gc/($gc+$gt)));
        print join ("\t", $locus, $gscore, $gc, $gt), "\n";
    }

    my $c = ($gc - $ec);
    my $t = ($gt - $et);

    if ($c + $t > 0) {
        my $score =  $c / ($c + $t);
        # print join ("\t",
        #             $locus,
        #             sprintf ("%g", $score),
        #             $c,
        #             $t,
        #         ), "\n";
    }
}



sub index_gff {
    my $gff_file = shift;

    my %exons = ();

    open my $GFF, '<', $gff_file or croak "Can't open $gff_file: $!";

  FEATURE:
    while (<$GFF>) {
        next if m/^\s*#/;
        chomp;

        my %locus = %{gff_read ($_)};
        
        my ($locus_id, $c, $t) = split /;/, $locus{attribute};
        $locus_id ||= q{.};
        
        if (!defined $c or !defined $t) {
            next FEATURE;
        }
        else {
            $locus_id =~ s/ID=//;
            $c =~ s/c=//;
            $t =~ s/t=//;
        }

        $exons{$locus_id}->[0] += $c;
        $exons{$locus_id}->[1] += $t;
    }
    return \%exons;
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
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/collapse_exons.pl $:
 $Id: collapse_exons.pl 249 2010-01-12 05:24:34Z psilva $:

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
