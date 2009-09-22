#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw /sum/;

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my $type = 'verbose';
my $frequencies;
my $id_regex;
my $reference;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'output|o=s'     => \$output,
    'type|t=s'       => \$type,
    'frequencies|f'  => \$frequencies,
    'id-regex|i=s'   => \$id_regex,
    'reference|r=s'  => \$reference,
    'verbose|v'      => sub { use diagnostics; },
    'quiet|q'        => sub { no warnings; },
    'help|h'         => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'       => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result
and ($type eq 'concise' xor $type eq 'verbose');


# redirect standard output to file if requested
if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}


# read in bowtie verbose file
my $counts   = undef;
my $previous = undef;

while (<>) {
    chomp;
 
    my $current = read_bowtie ($_);

    if ($frequencies) {

        $counts->{$current->{target}->[0]}{alternatives} += $current->{alternatives};
        $counts->{$current->{target}->[0]}{frequencies}++;
    }
    else {
        unless (defined $previous) {
            $previous = $current;
        }
        elsif ($current->{read_id} eq $previous->{read_id}) {
            push @{$previous->{strand}}, $current->{strand}->[0];
            push @{$previous->{target}}, $current->{target}->[0];
            push @{$previous->{coordinate}}, $current->{coordinate}->[0];
            push @{$previous->{snp}}, $current->{snp}->[0];
        }
        else {
            print_eland (
                $previous->{read_id},
                $previous->{sequence},
                $previous->{target},
                $previous->{coordinate},
                $previous->{strand},
                $previous->{snp},
            );
            $previous = $current;
        }
    }
}


count_reads ($reference, $counts, $id_regex)
if $frequencies;

print_eland (
    $previous->{read_id},
    $previous->{sequence},
    $previous->{target},
    $previous->{coordinate},
    $previous->{strand},
    $previous->{snp},
) if defined $previous;



### done


sub read_bowtie {
    my ($bowtie_line) = @_;

    my ($read_id, $strand, $target, $coordinate, $sequence, $qualities, $alternatives, $snp)
    = split /\t/, $bowtie_line;

    my @mm = (split /,/, $snp);

    @mm = () if $snp =~ m/^\s*$/;

    return {
        'read_id'      => $read_id,
        'strand'       => [$strand],
        'target'       => [$target],
        'coordinate'   => [$coordinate],
        'sequence'     => $sequence,
        'qualities'    => $qualities,
        'alternatives' => $alternatives,
        'snp'          => [scalar @mm],
    }
}


sub print_eland {
    my ($read_id, $sequence, $chromosomes_ref, $coordinates_ref, $strands_ref, $mismatches_ref)
    = @_;

    croak "Total number of chromosomes, coordinates, strands, and mismatches don't match"
    unless scalar @{$chromosomes_ref} == scalar @{$coordinates_ref} 
    and scalar @{$chromosomes_ref} == scalar @{$strands_ref}
    and scalar @{$chromosomes_ref} == scalar @{$mismatches_ref};

    my $target;
    map {
        $target
        .= $chromosomes_ref->[$_]
        .  q{:} . $coordinates_ref->[$_]
        .  ($strands_ref->[$_] eq q{+} ? q{F} : q{R})
        .  $mismatches_ref->[$_]
        .  ($_ < @{$chromosomes_ref} - 1 ? q{,} : q{})
    }
    ( 0 .. @{$chromosomes_ref} - 1 );


    my %rank = (
        0 => 0,
        1 => 0,
        2 => 0
    );

    for (@{$mismatches_ref}) {
        $rank{$_}++
    }

    map { die if $rank{$_} eq 'string' ; $rank{string} .= $rank{$_} . q{:} } sort keys %rank;
    chop $rank{string};

    print join ("\t",
                $read_id,
                $sequence,
                $rank{string},
                $target,
            ), "\n";
}


sub count_reads {
    my ($reference, $counts_ref, $id_regex) = @_;

    return unless $reference;

    # read in chromosome/model lengths
    my %reference = %{ index_fasta ($reference) };

    # sort and print target id, read frequency on mapped to target per kb, average alternative mappings
  TARGET:
    for my $target (sort keys %{$counts_ref}) {

        my ($id) = $target =~ m/$id_regex/;
print STDERR $id, "\n";
        unless (exists $reference{$target}) {
            carp "$target doesn't exist in $reference\n";
            next TARGET;
        }

        print join ("\t",
                    #$target,
                    $id,
                    sprintf ("%g", ($counts_ref->{$target}{frequencies} / length $reference{$target}) * 1000),
                    sprintf ("%g", ($counts_ref->{$target}{frequencies} ? $counts_ref->{$target}{alternatives} / $counts_ref->{$target}{frequencies} : 0)),
                ), "\n";
    }

}


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
            $fastaseq[$i] =~ s/[\r\n]//g;
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
 -f, --frequencies count read frequencies in the process
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
