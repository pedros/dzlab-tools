#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

# Check required command line parameters
unless (@ARGV > 0) {
    pod2usage ( -verbose => 1 );
}

my $reference;

# Grabs and parses command line options
my $result = GetOptions (
    'reference|r=s' => \$reference,
    'verbose|v' => sub { use diagnostics; },
    'quiet|q'   => sub { no warnings; },
    'help|h'    => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'  => sub { pod2usage ( -verbose => 2 ); }
);


my %reference = %{ index_fasta ($reference) };

while (<>) {
    chomp; # delete line feeds
    next if ($_ =~ m/^#.*$|^\s*$|\.;\.$/);

    my %eland = %{ parse_eland3_line($_) };

    for my $i (0 .. $eland{matches}) {
        print join ("\t",
                    $eland{chr},
                    q{.},
                    "$eland{line}:$eland{sequence}",
                    $eland{"coord$i"},
                    ($eland{"coord$i"} + length $eland{sequence}),
                    $eland{matches},
                    $eland{"strand$i"},
                    $eland{"mm$i"},
                    



    print ">$attribute $site{seqname}:$site{start}:$site{end}\n";
    print substr $reference{$site{seqname}}, ($site{start} - 1), ($length);
    print "\n";
}


# index_fasta reads a standard fasta file and returns a hash with chromosome names as keys as chromosome lengths as values
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


# Parses single eland3 record (output from seqmap)
# Returns reference to hash with following keys:
# 'line', 'sequence', 'matches', 'chr[0-#matches-1]', 'coord[0-#matches-1]'
sub parse_eland3_line {
    my $tmp_eland_line = shift;
    $tmp_eland_line =~ s/[\n\r]//g;
    my @eland_line = (split /\t/, $tmp_eland_line);

    my %hash = ();
    $hash{'line'}     = $eland_line[0];
    $hash{'sequence'} = $eland_line[1];

    if ($eland_line[2] =~ m/NM/) {
        $hash{'matches'} = 0;
    }

    elsif ($eland_line[2] =~ m/^[0-9]+$/) {
        $hash{'matches'} = $eland_line[2];
    }

    elsif ($eland_line[2] =~ m/:/) {
        my @all_matches = split /,/, $eland_line[3];
        $hash{'matches'} = scalar @all_matches;
    }

    if ($hash{'matches'} > 1) {
        my @all_reads = split /,/, $eland_line[3];

        for my $i (0..@all_reads - 1) {
            ($hash{"chr$i"}, $hash{"coord$i"}) = split /:/, $all_reads[$i];
            ($hash{'strand$i'},$hash{'mm$i'})  = hash{'coord$i'} =~ s/([A-Z])([0-9])$//i;
            $hash{'strand$i'} = ($hash{'strand$i'} =~ m/F/i) ? q{+} : q{-};
        }
    }
    elsif ($hash{'matches'} == 1) {
        ($hash{'chr0'}, $hash{'coord0'}) = split /:/, $eland_line[3];
        ($hash{'strand0'},$hash{'mm0'})  = hash{'coord0'} =~ s/([A-Z])([0-9])$//i;
        $hash{'strand0'} = ($hash{'strand0'} =~ m/F/i) ? q{+} : q{-};
    }
    return \%hash;
}



__END__


=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

=head1 REVISION

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
