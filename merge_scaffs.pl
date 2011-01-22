#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my $scaffolds;
my $prefix;
my $max_length = 30000000;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'fasta-scaffolds|f=s' => \$scaffolds,
    'max-length|l=i'      => \$max_length,
    'prefix|p=s'          => \$prefix,
    'output|o=s'          => \$output,
    'verbose|v' => sub { use diagnostics; },
    'quiet|q'   => sub { no warnings; },
    'help|h'    => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'  => sub { pod2usage ( -verbose => 2 ); }
);

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my $current_length = 0;
my @objects = ();

my @reference = @{ index_fasta ($scaffolds) };

my $last_offset;
my $last_header;

while (my $scaff = shift @reference) {

    if ($current_length + length $scaff->[1] < $max_length and @reference) {
        push @objects, $scaff;
        $current_length += length $scaff->[1];
    }
    else {
        my $header = "$objects[0]->[0]-$objects[$#objects]->[0]";

        if ($prefix) {
            $header =~ s/$prefix//g;
            $header = 'group_' . $header;
        }

        print ">$header\n";

        my $offset = 0;
        for my $i (0 .. $#objects) {
            print STDERR "$header\t$objects[$i]->[0]\t$offset\n";
            print $objects[$i]->[1], "\n";
            $offset += length $objects[$i]->[1];
        }

        $current_length = 0;
        @objects = ();

        push @objects, $scaff;
        $last_offset = $offset + length $scaff->[1];
        $last_header = $header;
        $current_length += length $scaff->[1];
    }
}

if (@objects) {
    print STDERR "$last_header\t$objects[0]->[0]\t$last_offset\n";
    print $objects[0]->[1], "\n"
}

exit 0;


# index_fasta reads a standard fasta file and returns a hash with chromosome names as keys as chromosome lengths as values
sub index_fasta {
    my $reference_file = shift;

    # holds name of chromosomes as keys and length of chromosomes in bp as values
    my @reference = ();

    return \@reference unless $reference_file;

    # reads in the reference genome file into @fastaseq
    open my $REF, '<', "$reference_file" or croak "Can't open file: $reference_file";
    my @fastaseq = <$REF>;
    close $REF;

    # find and store indices for each chromosome change and corresponding descriptions
    my ( @idx, @dsc ) = ();
    for my $i ( 0 .. @fastaseq - 1 ) {
        if ( $fastaseq[$i] =~ m/^>/ ) {
            $fastaseq[$i] =~ s/>//gi;
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
        push @reference, [$dsc[$j], $line];

    }
    return \@reference;
}



__END__


=head1 NAME

 merge_scaffs.pl -- Merge fasta scaffolds into pseudo chromosomes of given maximum length

=head1 SYNOPSIS

 ./merge_scaffs.pl -l 30000000 -f multiple_scaffs.fa -o pseudo-chr.fa
 ./merge_scaffs.pl -l 30000000 -f multiple_scaffs.fa -o pseudo-chr.fa --prefix scaffold_ 2> offsets

=head1 DESCRIPTION

=head1 OPTIONS

 -l, --max-length        Maximum length for each generated pseudo chromosome
 -f, --fasta-scaffolds   Fasta file with unknown scaffolds
 -p, --prefix            remove common prefix from scaffolds headers [scaffold_]
 -o, --output            Output pseudo chromosome fasta file name

=head1 REVISION

 $Rev: 249 $:
 $Author: psilva $:
 $Date: 2010-01-11 21:24:34 -0800 (Mon, 11 Jan 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/merge_scaffs.pl $:
 $Id: merge_scaffs.pl 249 2010-01-12 05:24:34Z psilva $:


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
