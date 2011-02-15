#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if $opt_help || !$opt_fasta_file;

use FindBin;
use lib "$FindBin::Bin/DZLab-Tools/lib";
use DZLab::Tools::Fasta;


if ($opt_output) {
    open my $USER_OUT, '>', $opt_output or croak "Can't open $opt_output for writing: $!";
    select $USER_OUT;
}

my %reference = %{ slurp_fasta ($opt_fasta_file) };

my $total_bp     = 0;
my $total_non_bp = 0;

if ($opt_list) {
    print $opt_fasta_file, ":\n";

    for (sort keys %reference) {

        my $tmp_length_bp     = length $reference{$_};
        my $tmp_length_non_bp = $reference{$_} =~ tr/ACGTacgt//c;

        print join ("\t",
                    $_,
                    $tmp_length_bp,
                    $tmp_length_bp - $tmp_length_non_bp,
                ), "\n";

        $total_bp += $tmp_length_bp;
        $total_non_bp += $tmp_length_non_bp;
    }
    print "Total size:\t$total_bp\t", $total_bp - $total_non_bp, "\n";

    exit 0;
}


if ($opt_seqid) {

    if (-e $opt_seqid) {
        
        open my $SEQID, '<', $opt_seqid or croak "Can't open $opt_seqid for reading";

        SEQID:
        while ($opt_seqid = <$SEQID>) {
            $opt_seqid =~ s/[\n\r]//g;

            my $sequence = _get_subseq (\%reference, $opt_seqid, $opt_rc);

            next SEQID unless $sequence;

            my ($start,$end) = %opt_range ? @opt_range{'start','end'} : (1, length $sequence);

            my $file_name = fileparse ($opt_fasta_file);

            # one time masking. keeping it here in case it's needed in the future
            # if ($opt_seqid =~ m/^scaffold_2$/) {
            #     print STDERR length $sequence, "\n";
            #     my $n = 'N' x (2350000 - 1650000);
            #     substr $sequence, (1650000 - 1), (2350000 - 1650000), $n;
            # }
            # elsif ($opt_seqid =~ m/^scaffold_175$/) {
            #     print STDERR length $sequence, "\n";
            #     my $n = 'N' x (56000 - 46000);
            #     substr $sequence, (46000 - 1), (56000 - 46000), $n;
            # }

            print ">$opt_seqid\n";
            print $sequence, "\n";
        }
    }
    else {

        my $sequence = _get_subseq (\%reference, $opt_seqid, $opt_rc, %opt_range);

        my ($start,$end) = %opt_range ? @opt_range{'start','end'} : (1, length $sequence);

        my $file_name = fileparse ($opt_fasta_file);

        print format_fasta("lcl|$file_name|$opt_seqid|$start-$end", $sequence), "\n";

        exit 0;
    }
}


if ($opt_split) {
    for my $chr (sort keys %reference) {
        open my $CHROUT, '>', "$opt_fasta_file-$chr" or croak "Can't write to $opt_fasta_file-$chr";

        print $CHROUT ">$chr\n";
        print $CHROUT $reference{$chr}, "\n";

        close $CHROUT;
    }

    exit 0;
}


sub _get_subseq {
    my ($seq_href, $seqid, $rc, $start, $end) = @_;

    $seqid =~ tr/A-Z/a-z/;

    unless (exists $seq_href->{$seqid}) {
        croak "Sequence ID $seqid does not exist";
    }

    my $sequence = $seq_href->{$seqid};

    if ($start and $end) {
        # b/c start, end are 1-based
        my $s0 = $start-1;
        my $e0 = $end-1;
        my $last = length($sequence) - 1;

        croak "Coordinates out of bounds" if ($s0 < 0 || $e0 > $last);

        my $subseq = substr ($sequence, $s0 , $e0 - $s0);
        $subseq =~ tr/ACGTacgt/TGCAtgca/ if $rc;
        return $subseq;
    }
    else {
        $sequence =~ tr/ACGTacgt/TGCAtgca/ if $rc;
        return $sequence;
    }
}


=head1 NAME

parse_fasta.pl - Retrieve sequence information from fasta files

=head1 SYNOPSIS

list all sequences in fasta file:

 parse_fasta.pl -l fasta.txt

print base pair 123 to 456 in sequence chr1:

 parse_fasta.pl -r 123 456 -i chr1 fasta.txt

retrieves, reverse compliments base pair 123 to 456 in sequence chr1
(NOTE: all coordinates are relative to the 5' end of the + strand):

 parse_fasta.pl --rc -r 123 456 -i chr1 fasta.txt

=head1 OPTIONS

=over 1

=item -l | --list        

print list of sequence ids and lengths

=item -i <sequence_name> | --seqid <sequence_name>  

sequence id from which to print sub sequence

=for Euclid
    sequence_name.type: string

=item -r <start> <end> | --range <start> <end>

start and end coordinates to print

=for Euclid
    start.type: int >=0
    end.type: int >=0

=item -s | --split       

split input fasta file into <basename> chromosomes

=item --rc          

take reverse compliment of range

=item -o <file> | --output <file>  

filename to write results to (defaults to STDOUT)

=for Euclid
    file.type: writeable

=item -v | --verbose     

output perl's diagnostic and warning messages

=item -q | --quiet       

supress perl's diagnostic and warning messages

=item -h | --help

show this message.

=item <fasta_file>

Name of fasta file.

=for Euclid
    fasta_file.type: readable

=back

=head1 REVISION

 Version 0.0.1

 $Rev: 249 $:
 $Author: psilva $:
 $Date: 2010-01-11 21:24:34 -0800 (Mon, 11 Jan 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/parse_fasta.pl $:
 $Id: parse_fasta.pl 249 2010-01-12 05:24:34Z psilva $:

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
