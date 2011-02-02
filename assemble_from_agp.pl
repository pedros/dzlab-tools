#!/usr/bin/env perl
# ___UNDOCUMENTED___

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my $output;
my $agp_file;
my $reference;

# Grabs and parses command line options
my $result = GetOptions (
    'agp-file|f=s'  => \$agp_file,
    'reference|r=s' => \$reference,
    'output|o=s'    => \$output,
    'verbose|v'     => sub { use diagnostics; },
    'quiet|q'       => sub { no warnings; },
    'help|h'        => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'      => sub { pod2usage ( -verbose => 2 ); }
);

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

my %reference = %{ index_fasta ($reference) };
my %buffer    = ();
my %objects   = ();

open my $AGP, '<', $agp_file or croak "Can't open $agp_file: $!";
while (<$AGP>) {

    my %agp = %{ read_agp ($_) };

    unless (exists $objects{$agp{object}}) {
        print "\n>$agp{object}\n";
        $objects{$agp{object}} = undef;
    }

    if ($agp{component_type} !~ m/N/i) {
        my $string = substr $reference{$agp{component_id}}, $agp{component_beg}, ($agp{component_end} - $agp{component_beg} + 1);
        $string = reverse_complement ($string) if $agp{orientation} !~ m/\+/;
#        $buffer{object} .= $string;
        print $string;
    }
    else {
#        $buffer{object} .= 'N' x $agp{gap_length};
        my $string = 'N' x $agp{gap_length};
        print $string;
    }

#    (print $buffer{object}) && delete $buffer{object}
#    if length $buffer{object} > 1000;
}
close $AGP or croak "Can't close $agp_file: $!";


sub read_agp {
    chomp;
    my ($object, $object_beg, $object_end, $part_number, $component_type,
        $component_id_or_gap_length, $component_beg_or_gap_type, $component_end_or_linkage, $orientation_or_empty)
    = split /\t/;

    $object =~ tr/A-Z/a-z/;
    $component_id_or_gap_length =~ tr/A-Z/a-z/;
    $component_beg_or_gap_type--;
    $component_end_or_linkage--;

    my %agp = (
	'object'         => $object,
	'object_beg'     => $object_beg,
	'object_end'     => $object_end,
	'part_number'    => $part_number,
	'component_type' => $component_type,
    );

    if ($agp{component_type} =~ m/N/) {
        $agp{gap_length} = $component_id_or_gap_length;
        $agp{gap_type}   = $component_beg_or_gap_type;
        $agp{linkage}    = $component_end_or_linkage;
    }
    else {
        $agp{component_id}  = $component_id_or_gap_length;
        $agp{component_beg} = $component_beg_or_gap_type;
        $agp{component_end} = $component_end_or_linkage;
        $agp{orientation}   = $orientation_or_empty;
    }

    return \%agp;
}


sub print_fasta {
    my $sequence = shift;
    my $line_max = shift || 80;

    $sequence =~ s/\r\n//;

    my $i;
    for ($i = 0; $i < (length $sequence) - $line_max - 1; $i += $line_max) {
        print substr $sequence, $i, $line_max - 1;
        print "\n";
    }
    print substr $sequence, $i;
    print "\n";

}


sub reverse_complement {
    my $sequence = shift;
    $sequence =~ tr/ACGTacgt/TGCAtgca/;
    return reverse $sequence;
}


sub index_fasta {
    my $reference_file = shift;

    # holds name of chromosomes as keys and length of chromosomes in bp as values
    my %reference = ();

    return \%reference unless $reference_file;

    # reads in the reference genome file into @fastaseq
    open my $REF, '<', "$reference_file" or croak "Can't open file: $reference_file";
    my @fastaseq = <$REF>;
    close $REF;

    # find and store indices for each chromosome change and corresponding descriptions
    my ( @idx, @dsc ) = ();
    for my $i ( 0 .. @fastaseq - 1 ) {
        if ( $fastaseq[$i] =~ m/^>/ ) {
            $fastaseq[$i] =~ s/>//g;
            $fastaseq[$i] = ( split /\s/, "$fastaseq[$i]" )[0];
            $fastaseq[$i] =~ tr/A-Z/a-z/;
            push @idx, $i;
            push @dsc, $fastaseq[$i];
        }
    }

    # tries to find each chromosome's centrometer center coordinate
    for my $j ( 0 .. @idx - 1 ) {
        my $line;
        if ( $j == scalar @idx - 1 ) {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. @fastaseq - 1]);
        }
        else {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. $idx[$j + 1] - 1]);
        }
        $line =~ s/[\n\r]//g;
        $reference{$dsc[$j]} = $line;
    }
    return \%reference;
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
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/assemble_from_agp.pl $:
 $Id: assemble_from_agp.pl 249 2010-01-12 05:24:34Z psilva $:

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
