#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

my $output;
my $feature;
my $sequence;

# Grabs and parses command line options
my $result = GetOptions (
    'feature|f=s'  => \$feature,
    'sequence|s=s' => \$sequence,
    'output|o=s'   => \$output,
    'verbose|v'  => sub { use diagnostics; },
    'quiet|q'    => sub { no warnings; },
    'help|h'     => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'   => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result and -e $ARGV[0]
or (
    ($feature and !$sequence) or
    (!$feature and $sequence)
);

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

my ($name, $path, $suffix) = fileparse ($ARGV[0], qr/\.[^.]*/);

my %file_handle;
while (my $gff_line = <>) {
    next if $gff_line =~ m/^\s*$|^\s*#/;

    my $current = (split /\t/, $gff_line)[($sequence ? 0 : 2)];

    next unless $feature =~ m/all/i
    or $current =~ m/$feature/i;

    my $out = $path . $name . "-$current" . $suffix;

    open $file_handle{$current}, '>', $out
    or croak "Can't open $out for writing: $!"
    unless exists $file_handle{$current};

    print {$file_handle{$current}} $gff_line;
}
close $file_handle{$_} for keys %file_handle;

exit 0;


__END__


=head1 NAME

 split_gff - Split GFF files by sequence ID or feature

=head1 SYNOPSIS

 split_gff.pl -f exon all_features.gff  # filter by exon
 split_gff.pl -s chr1 all_sequences.gff # filter by chromosome
 split_gff.pl -f all al_features.gff    # create multiple files, one per feature

=head1 DESCRIPTION

=head1 OPTIONS

 split_gff.pl [OPTION]... [FILE]...

 -f, --feature     feature used to filter GFF file by ('all' generates one file per feature)
 -s, --sequence    sequence ID used to filter GFF file by ('all' generates one file per sequence ID)
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
