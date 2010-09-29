#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

my $DATA_HANDLE = 'ARGV';
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'output|o=s'  => \$output,
    'verbose|v'   => sub { use diagnostics; },
    'quiet|q'     => sub { no warnings; },
    'help|h'      => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'    => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $result;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my @buffer = ();
while (<$DATA_HANDLE>) {
    next if /\s*#/;
    my @fields = split /\t/;

    if ( !@buffer || ($fields[3] - 1 == $buffer[-1]->[4] and $fields[5] == $buffer[-1]->[5]) ) {
        push @buffer, [@fields[0,1,2], 1, @fields[4,5,6,7,8]] if !@buffer;
        push @buffer, \@fields;
    }
    else {
        print join ("\t", @{$buffer[0]}[0,1,2,3], $buffer[-1]->[4], @{$buffer[0]}[5,6,7,8]);
        print join ("\t", @{$buffer[0]}[0,1,2], $buffer[-1]->[4] + 1, $fields[4] - 1, 0, @{$buffer[0]}[6,7,8]);
        @buffer = (\@fields);
    }
}

print join ("\t", @{$buffer[0]}[0,1,2,3], $buffer[-1]->[4], @{$buffer[0]}[5,6,7,8])
if @buffer;


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
