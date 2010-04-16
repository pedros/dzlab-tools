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
my $paired;
my $eland;
my $id_regex;
my $reference;
my $output;
my $unmatched;
my @splice;

# Grabs and parses command line options
my $result = GetOptions (
    'output|o=s'     => \$output,
    'recover|u=s'    => \$unmatched,
    'splice|s=i{2}'  => \@splice,
    'type|t=s'       => \$type,
    'frequencies|f'  => \$frequencies,
    'paired|p'       => \$paired,
    'eland|e'        => \$eland,
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
and ($type eq 'concise' xor $type eq 'verbose')
and (   ($frequencies xor $paired)
     or (!$frequencies and !$paired));

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
    s/[\r\n]//;

    my $current = read_bowtie ($_);

    $current->{snps}{$current->{snp}->[0]}++;

    if ($frequencies) {
        $counts->{$current->{target}->[0]}{alternatives}
        += $current->{alternatives};
        $counts->{$current->{target}->[0]}{frequencies}++;
    }
    elsif ($paired) {

        my $next = <>; chomp $next; $next =~ s/[\r\n]//;
           $next = read_bowtie ($next);

        $next->{snps}{$next->{snp}->[0]}++;

        print_gff ($current, $next);
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
            $previous->{snps}{$current->{snp}->[0]}++;
        }
        else {

            catch_up ($previous, $unmatched, @splice)
            if $unmatched;

            print_eland ($previous);
            $previous = $current;
        }
    }
}

count_reads ($reference, $counts, $id_regex) if $frequencies;

catch_up ($previous, $unmatched, @splice) if defined $previous and $unmatched;

print_eland ($previous) if defined $previous;

# for when last read in bowtie file is *not* last read in fasta file
catch_up ($previous, $unmatched, @splice) if defined $previous and $unmatched;



{
    my %file_handles;
    my $file_handle;

    sub catch_up {
        my ($current, $unmatched, @splice) = @_;

        $file_handle = $file_handles{$unmatched};

        unless (defined $file_handle) {
            open $file_handle, '<', $unmatched
            or croak "Can't open $unmatched: $!";
            $file_handles{$unmatched} = $file_handle;
        }

      FASTA_HEADER:
        while (    defined (my $header   = <$file_handle>)
               and defined (my $sequence = <$file_handle>)) {

            chomp $header; chomp $sequence;

            $header =~ s/^([>@])//;
            if ( q{@} eq $1 ) {
                <$file_handle>;
                <$file_handle>;
            }
            elsif ( q{>} ne $1 ) {
                croak "Can't figure out whether this file is straight fasta or fastq"
            }

            my %unmatched = (header => $header, sequence => $sequence);

            # if potentially unmatched read is not current bowtie read
            if ($unmatched{header} !~ m/$current->{read_id}/) {

                $unmatched{sequence}
                = substr $unmatched{sequence}, ($splice[0] - 1), ($splice[1] - $splice[0] + 1)
                if @splice;

                print join ("\t",
                            $unmatched{header},
                            $unmatched{sequence},
                            'NM',
                            "\n"
                        );
            }
            else {
                $current->{sequence} = $sequence;

                $current->{sequence}
                = substr $current->{sequence}, ($splice[0] - 1), ($splice[1] - $splice[0] + 1)
                if @splice;

                last FASTA_HEADER;
            }
        }
    }
}

sub print_gff {
    my ($current, $next) = @_;

    print join ("\t",
                $current->{target}->[0],
                'bowtie',
                'frag',
                $current->{coordinate}->[0],
                $next->{coordinate}->[0] + length ($next->{sequence}) - 1,
                q{.},
                q{+},
                q{.},
                "alt=$current->{alternatives}"
            ), "\n";
}


sub print_eland {
    my ($previous) = @_;

    my ($read_id, $sequence, $chromosomes_ref, $coordinates_ref, $strands_ref, $mismatches_ref, $mismatches_total)
    = ($previous->{read_id}, $previous->{sequence}, $previous->{target}, $previous->{coordinate}, $previous->{strand}, $previous->{snp}, $previous->{snps});

    croak "Total number of chromosomes, coordinates, strands, and mismatches don't match"
    unless scalar @{$chromosomes_ref} == scalar @{$coordinates_ref} 
    and    scalar @{$chromosomes_ref} == scalar @{$strands_ref}
    and    scalar @{$chromosomes_ref} == scalar @{$mismatches_ref};

    my $mismatches
    = join q{:}, map { $mismatches_total->{$_} || 0 } (0 .. 2);

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

    print join ("\t",
                $read_id,
                $sequence,
                $mismatches,
                $target,
            ), "\n";
}


sub read_bowtie {
    my ($bowtie_line) = @_;

    my ($read_id, $strand, $target, $coordinate, $sequence, $qualities, $alternatives, $snp)
    = split /\t/, $bowtie_line;

    my @mm = split /,/, $snp;

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
 -u, --recover     given original alignment fasta file, recovers unmatched reads (which bowtie does not output)
 -p, --paired      convert bowtie's paired ends output to gff with concatenated library ends per region
 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.2

 $Rev$:
 $Author$:
 $Date$:
 $HeadURL$:
 $Id$:

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
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
