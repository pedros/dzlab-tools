#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

my $offsets;
my $new_feature;
my $upstream   = 0;
my $downstream = 0;
my $skip_absent;
my $trim;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'offsets|f=s'     => \$offsets,
    'trim|t'          => \$trim,
    'new-feature|n=s' => \$new_feature,
    'upstream|u=f'    => \$upstream,
    'downstream|d=f'  => \$downstream,
    'skip-absent|s'   => \$skip_absent,
    'output|o=s'      => \$output,
    'verbose|v'       => sub { use diagnostics; },
    'quiet|q'         => sub { no warnings; },
    'help|h'          => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'        => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $offsets xor ($upstream or $downstream);


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

LOCUS:
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
    elsif ($trim) {
        my ($start, $end) = @gff[3,4];
        $gff[3] = $start + $downstream;
        $gff[4] = $end   - $upstream;

	$gff[3] = 0 if $gff[3] < 0;
        next LOCUS unless $gff[3] <= $gff[4];
    }
    else {
        my ($start, $end) = @gff[3,4];

        if ($gff[6] eq q{-}) {
            $gff[3] = $end - $downstream;
            $gff[4] = $end + $upstream;
        }
        else {
            $gff[3] = $start - $upstream;
            $gff[4] = $start + $downstream;
        }
        
	$gff[3] = 0 if $gff[3] < 0;

        next LOCUS unless $gff[3] <= $gff[4];
    }

    $gff[2] = $new_feature if $new_feature;

    print join "\t", @gff;
}


sub gff_read {
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute)
    = split /\t/, shift;

    return {
	'seqname'  => $seqname,
	'source'   => $source,
	'feature'  => $feature,
	'start'    => $start,
	'end'      => $end,
	'score'    => $score,
	'strand'   => $strand,
	'frame'    => $strand,
	'attribute'=> $attribute
    };
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
 -t, --trim        trim -u from upstream and -d from downstream
 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: 295 $:
 $Author: psilva $:
 $Date: 2010-04-05 16:39:12 -0700 (Mon, 05 Apr 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/offset_coordinates.pl $:
 $Id: offset_coordinates.pl 295 2010-04-05 23:39:12Z psilva $:

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
