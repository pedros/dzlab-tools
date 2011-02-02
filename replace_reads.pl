#!/usr/bin/env perl
# ___UNDOCUMENTED___

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

my $output;
my $fasta;
my $readsize;
my @splice;

# Grabs and parses command line options
my $result = GetOptions (
    'fasta|f=s'           => \$fasta,
    'readsize|r=i'        => \$readsize,
    'splice|s=i{,}'       => \@splice,
    'output|o=s'          => \$output,
    'verbose|v'           => sub { use diagnostics; },
    'quiet|q'             => sub { no warnings; },
    'help|h'              => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'            => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $result and @ARGV;

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

my $left_offset  = $splice[0] - 1;
#my $right_offset = $readsize - $splice[1];

open my $FASTA, '<', $fasta or croak "Can't open $fasta: $!";
while (<>) {
    chomp;
    my @fields = split /\t/;

    my $header   = <$FASTA>;
    my $sequence = <$FASTA>;
    $sequence =~ s/[\r\n]//g;

    croak "Read IDs don't match:\n$fields[0]\n$header"
    unless $header =~ m/${fields[0]}/;

    if ( (length $sequence) != (length $fields[1]) ) {
        $sequence    = substr $sequence, ($splice[0] - 1), (length $fields[1]);
        $left_offset = 0;
    }

    $fields[1] = $sequence;

    $fields[3] =~ s/(\w+:)(\d+)?(\w+)/$1.($2 - $left_offset).$3/eg
    unless $fields[2] =~ 'NM' or $left_offset == 0;

    print join ("\t", @fields), "\n";
}
close $FASTA or warn "Can't close $fasta: $!";


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
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/replace_reads.pl $:
 $Id: replace_reads.pl 249 2010-01-12 05:24:34Z psilva $:

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
