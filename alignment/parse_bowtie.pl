#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my $type = 'verbose';
my $id_regex;
my $reference;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'output|o=s'    => \$output,
    'type|t=s'      => \$type,
    'id-regex|i=s'  => \$id_regex,
    'reference|r=s' => \$reference,
    'verbose|v'     => sub { use diagnostics; },
    'quiet|q'       => sub { no warnings; },
    'help|h'        => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'      => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result and ($type eq 'concise' or $type eq 'verbose') and $reference;

# redirect standard output to file if requested
if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

# read in bowtie verbose file
my %counts = ();
while (<>) {
    chomp;
    
    my ($read_id, $strand, $target, $coordinate, $sequence, undef, $alternatives, $snp)
    = split /\t/;

    $target = (split /\s/, $target)[0];

    $counts{$target}{alternatives} += $alternatives;
    $counts{$target}{frequencies}++;
}

# read in chromosome/model lengths
my %reference = %{ index_fasta ($reference) };

# sort and print target id, read frequency on mapped to target per kb, average alternative mappings
TARGET:
for my $target (sort keys %counts) {

    my ($id) = $target =~ m/$id_regex/;

    print STDERR "$target doesn't exist in $reference\n" && next TARGET
    unless exists $reference{$target};

    print join ("\t",
                #$target,
                $id,
                sprintf ("%g", ($counts{$target}{frequencies} / length $reference{$target}) * 1000),
                sprintf ("%g", ($counts{$target}{frequencies} ? $counts{$target}{alternatives} / $counts{$target}{frequencies} : 0)),
            ), "\n";
}

### done



sub index_fasta {
    my $reference_file = shift;

    my %reference = ();

    return \%reference unless $reference_file;

    # reads in the reference genome file into @fastaseq
    open my $REF, '<', "$reference_file" or croak "Can't open $reference for reading: $!";
    my @fastaseq = <$REF>;
    close $REF;

    # find and store indices for each chromosome change and corresponding descriptions
    my ( @idx, @dsc ) = ();
    for my $i ( 0 .. @fastaseq - 1 ) {
        if ( $fastaseq[$i] =~ m/^>/ ) {
            $fastaseq[$i] =~ s/>//g;
            $fastaseq[$i] = ( split /\s/, "$fastaseq[$i]" )[0];
#            $fastaseq[$i] =~ tr/A-Z/a-z/;
            push @idx, $i;
            push @dsc, $fastaseq[$i];
        }
    }

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

 parse_bowtie.pl - Process bowtie alignment file into simple fields

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

 parse_bowtie.pl [OPTION]... [FILE]...

 -t, --type        type of bowtie output file (verbose or concise -- only verbose supported for now)
 -i, --id-regex    perl-type regular expression to identify feature id (ie. gene) in fasta alignment header (must include capturing parenthesis)
 -r, --reference   genome/cDNA models file in fasta format (for calculating relative frequency scores, etc.)
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
