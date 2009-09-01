#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

my $offsets;
my $new_feature;
my $upstream;
my $downstream;
my $skip_absent;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'offsets|f=s'     => \$offsets,
    'new-feature|n=s' => \$new_feature,
    'upstream|u=i'    => \$upstream,
    'downstream|d=i'  => \$downstream,
    'skip-absent|s'   => \$skip_absent,
    'output|o=s'      => \$output,
    'verbose|v'       => sub { use diagnostics; },
    'quiet|q'         => sub { no warnings; },
    'help|h'          => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'        => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $offsets xor ($upstream and $downstream);


if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

my %offsets = ();

if ($offsets) {
    open my $OFFSETS, '<', $offsets or croak "Can't open $offsets: $!";
    while (<$OFFSETS>) {
	chomp;
	my ($group, $scaffold, $offset)
	    = split /\t/;
	
	$scaffold =~ tr/A-Z/a-z/;
	
	$offsets{$scaffold}{group} =  $group;
	$offsets{$scaffold}{offset} = $offset;
    }
    close $OFFSETS or croak "Can't close $offsets: $!";
}

while (<>) {
    next if ($_ =~ m/^#.*$|^\s*$/);
    my @gff = split /\t/, $_;

    $gff[0] =~ tr/A-Z/a-z/;

    if ($offsets) {
	unless (exists $offsets{$gff[0]}) {
	    warn "Scaffold $gff[0] doesn't exist in offsets file $offsets";
            next if $skip_absent;
	}
	else {
	    $gff[3] += $offsets{$gff[0]}{offset};
	    $gff[4] += $offsets{$gff[0]}{offset};
	    $gff[0] =  $offsets{$gff[0]}{group};
	}
    }
    else {
	my $five_prime;
	($gff[5] eq q{+}) ? $five_prime = $gff[3] : $five_prime = $gff[4];

	$gff[3] = $five_prime - $upstream;
	$gff[4] = $five_prime + $downstream;

	$gff[3] = 0 if $gff[3] < 0;
    }

    $gff[2] = $new_feature if $new_feature;

    print join "\t", @gff;
}


__END__


=head1 NAME

 offset_coordinates.pl - Change coordinates in a gff annotation file

=head1 SYNOPSIS

 # with a list of pre-computed offsets per scaffold/sequence
 offset_coordinates --offsets file-with-offset-list.dat -o new-annotation-file.gff old-annotation-file

 # to convert gene list to (predicted) promoter list
 offset_coordinates --upstream 500 --downstream 200 --new-feature promoter -o new-annotation-file.gff old-annotation-file

=head1 DESCRIPTION

=head1 OPTIONS

 offset_coordinates.pl [OPTION]... [FILE]...

 -f, --offsets     filename of offset file (group	scaffold	offset)
 -n, --new-feature new feature to subsitute in field 3
 -u, --upstream    offset upstream of 5' (incompatible with --offsets)
 -d, --downstream  offset downstream of 5' (incompatible with --offsets)
 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

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
