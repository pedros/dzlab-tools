#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw /max/;

my $output;
my $html = 0;

# Grabs and parses command line options
my $result = GetOptions (
    'html|f'              => \$html,
    'output|o=s'          => \$output,
    'verbose|v'           => sub { use diagnostics; },
    'quiet|q'             => sub { no warnings; },
    'help|h'              => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'            => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result;

my $a_fasta   = index_fasta ($ARGV[0]);
my $b_fasta   = index_fasta ($ARGV[1]);
my @a_headers = sort keys %$a_fasta;
my @b_headers = sort keys %$b_fasta;

open my $USER_OUT_A, '>', "$ARGV[0].snp.fa" or croak "Can't open $ARGV[1].snp.fa for writing: $!";
open my $USER_OUT_B, '>', "$ARGV[1].snp.fa" or croak "Can't open $ARGV[1].snp.fa for writing: $!";

for my $a_seqid (@a_headers) {

    my $b_seqid = shift @b_headers;

    my $a_string  = [split //, $a_fasta->{$a_seqid}];
    my $b_string  = [split //, $b_fasta->{$b_seqid}];

    my $positions = vector_difference ($a_string, $b_string);
    my $totals    = 0;

    my $sequence = q{};

    no warnings; # array indices might go out of bounds for different-sized strings
    $sequence .= $positions->[$_] ? $a_string->[$_] : uc ($a_string->[$_])
    for (0 .. @$positions - 1);
    select $USER_OUT_A; print &string_to_fasta ($sequence, $a_seqid);
    

    $sequence  = q{};
    $sequence .= $positions->[$_] ? $b_string->[$_] : uc ($b_string->[$_])
    for (0 .. @$positions - 1);
    select $USER_OUT_B; print &string_to_fasta ($sequence, $b_seqid);

    $totals += $positions->[$_] ? 0 : 1
    for (0 .. @$positions - 1);
    
    print STDERR "$totals mismatch[es] found on $a_seqid versus $b_seqid\n";
}

close $USER_OUT_A;
close $USER_OUT_B;


sub file_to_string {
    my $file = shift;

    open my $INPUTFILE, '<', $file or croak "Can't open $file";
    my @string = <$INPUTFILE>;
    close $INPUTFILE;

    return [ map { next if m/^\s*#/; chomp; m/[^\r\n]/; split // } @string ];
}


sub vector_difference {
    my ($a_ref, $b_ref) = @_;
    warn "strings have different sizes: ", scalar @$a_ref, q{, }, scalar @${b_ref}
    if @$a_ref != @$b_ref;

    my @positions = ();

    no warnings;
    $positions[$_] = $a_ref->[$_] eq $b_ref->[$_]
    for (0 .. max (scalar @$a_ref, scalar @$b_ref) - 1);

    return \@positions;
}


sub index_fasta {
    my $reference_file = shift;

    my %reference = ();

    return \%reference unless $reference_file;

    # reads in the reference genome file into @fastaseq
    open my $REF, '<', "$reference_file" or croak "Can't open $reference_file for reading: $!";
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
        $reference{$dsc[$j]} = lc $line;
    }
    return \%reference;
}


sub string_to_fasta {
    my ($sequence, $header, $html) = @_;
    return unless $sequence;
    my $fasta_length   = 80;
    my $fasta_sequence = q{};
    
    $fasta_sequence .= ">$header\n" if $header;
    
    my $buffer_size = 0;
    for (ref $sequence eq 'ARRAY' ? @$sequence : split //, $sequence) {

        if ($buffer_size == $fasta_length - 1) {
            $fasta_sequence .= "$_\n";
            $buffer_size = 0;
        }
        else {
            $fasta_sequence .= $_;
            $buffer_size++
        }
    }
    return $fasta_sequence, "\n";
}


__END__


=head1 NAME

 simple_diff.pl - Compute Hamming distance between two strings or two files and point out where the difference is

=head1 SYNOPSIS

 # compare two strings
 simple_diff.pl string_a string_b

 # compare two files
 simple_diff.pl file_a file_b

=head1 DESCRIPTION

=head1 OPTIONS

 simple_diff.pl [OPTION]... [FILE A]    [FILE B]
 simple_diff.pl [OPTION]... [STRING A] [STRING B]

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
