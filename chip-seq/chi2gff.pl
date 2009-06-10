#!/usr/bin/perl

use strict;
use warnings;

my $infile  = $ARGV[0];
my $logfile = $infile;
$logfile =~ s/peaks/STDOUT/;

open my $LOG, '<', $logfile or die "Can't open $logfile";
my @log = <$LOG>;
close $LOG or die "Can't close $logfile";

@log = grep {
    m/alpha.*=/
} @log;

my $pattern = '[+-]?\.*(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?';

my @feature
= split /,/, $log[0];

($feature[0]) = $feature[0] =~ m/($pattern)/;
($feature[1]) = $feature[1] =~ m/($pattern)/;

$feature[2]
= 'a=' . sprintf("%g", $feature[0]) . ';t=' . sprintf("%g", $feature[1]);

while (<>) {
    my @fields = split /\t/, $_;
    print join ("\t",
                $fields[1],
                q{.},
                $feature[2],
                $fields[2],
                $fields[3],
                $fields[4],
                q{.},
                q{.},
                $fields[5],
                "\n"
            );
}


__END__


=head1 NAME

 chi2gff.pl -- Converts ChIPotle output files to gff format

=head1 SYNOPSIS

 ./chi2gff.pl chipotle-out.peaks > chipotle-out.gff

=head1 DESCRIPTION

 Assumes chipotle-out.STDOUT is in same directory as input file, which it uses to extract the parameters given to chipotle
 to use in the feature field.

=head1 OPTIONS

=head1 REVISION

 Version 0.0.1

 $Rev$:
 $Author$:
 $Date$:
 $HeadURL$:
 $Id$:

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
