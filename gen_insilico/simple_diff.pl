#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw /max/;

my $output;
# Grabs and parses command line options
my $result = GetOptions (
    'output|o=s'          => \$output,
    'verbose|v'           => sub { use diagnostics; },
    'quiet|q'             => sub { no warnings; },
    'help|h'              => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'            => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}
    
my $a_string = (-e $ARGV[0] ? file_to_string ($ARGV[0]) : [split //, $ARGV[0]]);
my $b_string = (-e $ARGV[1] ? file_to_string ($ARGV[1]) : [split //, $ARGV[1]]);
my $positions = vector_difference ($a_string, $b_string);
my $totals = 0;

no warnings; # array indices might go out of bounds for different-sized strings
print $a_string->[$_]
for (0 .. @$positions - 1);

print "\n";
print $positions->[$_] ? q{ } : q{|}
for (0 .. @$positions - 1);

print "\n";
print $b_string->[$_]
for (0 .. @$positions - 1);

$totals += $positions->[$_] ? 0 : 1
for (0 .. @$positions - 1);

print "\n";
print STDERR "$totals mismatch[es] found\n";


sub file_to_string {
    my $file = shift;

    open my $INPUTFILE, '<', $file or croak "Can't open $file";
    my @string = <$INPUTFILE>;
    close $INPUTFILE;

    return [ map { next if m/^\s*#/; chomp; m/[^\r\n]/; split // } @string ];
}


sub vector_difference {
    my ($a_ref, $b_ref) = @_;
    warn "strings have different sizes: ", scalar @$a_ref, q{, }, scalar @${b_ref}
    if @$a_ref != @$b_ref;

    my @positions = ();

    no warnings;
    $positions[$_] = $a_ref->[$_] eq $b_ref->[$_]
    for (0 .. max (scalar @$a_ref, scalar @$b_ref) - 1);

    return \@positions;
}

__END__


=head1 NAME

 simple_diff.pl - Compute Hamming distance between two strings or two files and point out where the difference is

=head1 SYNOPSIS

 # compare two strings
 simple_diff.pl string_a string_b

 # compare two files
 simple_diff.pl file_a file_b

=head1 DESCRIPTION

=head1 OPTIONS

 simple_diff.pl [OPTION]... [FILE A]    [FILE B]
 simple_diff.pl [OPTION]... [STRING A] [STRING B]

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
