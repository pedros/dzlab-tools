#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use feature 'say';

use FindBin;
use lib "$FindBin::Bin/DZLab-Tools/lib";
use DZLab::Tools::GFF qw/gff_to_string gff_make_iterator/;

my $feature;
my $sequence;

# Grabs and parses command line options
my $result = GetOptions (
    'feature|f=s'  => \$feature,
    'sequence|s=s' => \$sequence,
    'verbose|v'  => sub { use diagnostics; },
    'quiet|q'    => sub { no warnings; },
    'help|h'     => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'   => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV && $result && -e $ARGV[0] && ($feature xor $sequence);


my ($name, $path, $suffix) = fileparse ($ARGV[0], qr/\.[^.]*/);

my %file_handle;
my $iter = gff_make_iterator(file => $ARGV[0]);
RECORD:
while (defined(my $gff = $iter->())){
    next unless ref $gff eq 'HASH';

    my $current = $sequence ? $gff->{sequence} : $gff->{feature};

    if (($feature && ($feature eq 'all' || $feature eq $gff->{feature}))
        ||
        ($sequence && ($sequence eq 'all' || $sequence eq $gff->{sequence})))
    {
        my $out = $path . $name . "-$current" . $suffix;

        if (! exists $file_handle{$current}){
            open $file_handle{$current}, '>', $out
                or croak "Can't open $out for writing: $!";
        }

        say {$file_handle{$current}} gff_to_string $gff;
    }
}

__END__

=head1 NAME

 split_gff - Split GFF files by sequence ID or feature

=head1 SYNOPSIS

 split_gff.pl -f exon all_features.gff  # filter by exon, create file all_features-exon.gff
 split_gff.pl -s chr1 all_sequences.gff # filter by chromosome, create file all_sequences-chr1.gff
 split_gff.pl -f all all_features.gff   # create multiple files, one per feature

=head1 DESCRIPTION

=head1 OPTIONS

 split_gff.pl [OPTION]... [FILE]...

 -f, --feature     feature used to filter GFF file by ('all' generates one file per feature)
 -s, --sequence    sequence ID used to filter GFF file by ('all' generates one file per sequence ID)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: 249 $:
 $Author: psilva $:
 $Date: 2010-01-11 21:24:34 -0800 (Mon, 11 Jan 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/split_gff.pl $:
 $Id: split_gff.pl 249 2010-01-12 05:24:34Z psilva $:

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
