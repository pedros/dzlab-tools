#!/usr/bin/env perl

## TODO: keep mm calls difference for both reads as a measure of quality

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use autodie;

my $a_name;
my $b_name;

# Grabs and parses command line options
my $result = GetOptions (
    'a-ecotype|a=s'    => \$a_name,
    'b-ecotype|b=s'    => \$b_name,
    'verbose|v'   => sub { use diagnostics; },
    'quiet|q'     => sub { no warnings; },
    'help|h'      => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'    => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result and $a_name and $b_name;

open my $A_IN, '<', $ARGV[0];
open my $B_IN, '<', $ARGV[1];

unlink $ARGV[0], $ARGV[1];

open my $A_NAME, '>', $ARGV[0];
open my $B_NAME, '>', $ARGV[1];

open my $A_original, '>', $ARGV[0] . ".orig";
open my $B_original, '>', $ARGV[1] . ".orig";

CMP:
while (    defined (my $a_record = <$A_IN>)
        and defined (my $b_record = <$B_IN>)) {
    print $A_original $a_record;
    print $B_original $b_record;

    my ($a_id, $a_mm) = (split /\t/, $a_record)[0,2];
    my ($b_id, $b_mm) = (split /\t/, $b_record)[0,2];

    $a_mm = get_score ($a_mm);
    $b_mm = get_score ($b_mm);

    if (! defined $a_mm and ! defined $b_mm) {
        next CMP; # no matches at all
    }
    elsif (defined $a_mm and defined $b_mm) {
        next CMP if $a_mm == $b_mm;
        print $A_NAME $a_record if $a_mm < $b_mm;
        print $B_NAME $b_record if $a_mm > $b_mm;
    }
    elsif (! defined $a_mm) {
        print $B_NAME $b_record;
    }	
    elsif (! defined $b_mm) {
        print $A_NAME $a_record;
    }
    else {croak "Impossible situation:\n$a_record\n$b_record"}
}

close $A_NAME; close $B_NAME;
close $A_IN;   close $B_IN;

sub get_score {
    my ($mm) = @_;

    return if 'NM' eq $mm;
    my @mm = split /:/, $mm;

    for my $i (0 .. @mm - 1) {
        return $i if 1 == $mm[$i];
    }
}


__END__


=head1 NAME

 split_on_mismatches.pl - Filter sorted bowtie inputs into two files based on mismatch counts

=head1 SYNOPSIS

 # overwrites input files with correct imprinting alignments
 perl split_on_mismatches.pl -a Col -b Ler CxL_En_WT_Col.bowtie CxL_En_WT_Ler.bowtie

=head1 DESCRIPTION

=head1 OPTIONS

 split_on_mismatches.pl [OPTION]... [FILE]...

 -a, --a-ecotype   ecotype of first input file
 -b, --b-ecotype   ecotype of second input file
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
