#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

my $DATA_HANDLE = 'ARGV';
my $output;
my @fasta;
my @motifs;

# Grabs and parses command line options
my $result = GetOptions (
    'fasta|f=s{,}'  => \@fasta,
    'motifs|m=s{,}' => \@motifs,
    'output|o=s'  => \$output,
    'verbose'   => sub { use diagnostics; },
    'quiet'     => sub { no warnings; },
    'help'      => sub { pod2usage ( -verbose => 1 ); },
    'manual'    => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $result and @fasta and @motifs;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my %frequencies = ();

for my $motif (@motifs) {

    for my $fasta (@fasta) {

        my $fasta_iterator = make_fasta_iterator ($fasta);

        while (my $locus = $fasta_iterator->()) {
            
            # this should probably be inlined
            $frequencies{$motif}{$fasta}{$locus->{header}}
                = count_sites ($locus->{sequence}, $motif);

            print join ("\t",
                        $motif,
                        $fasta,
                        $locus->{header},
                        $frequencies{$motif}{$fasta}{$locus->{header}}
                    ), "\n";
        }

    }
}



# sub report_frequencies {
#     my ($frequencies) = @_;

#     for my $motif (sort keys %$frequencies) {
#         for 

# }

sub count_sites {
    my ($sequence, $type) = @_;

    my @count = $sequence =~ m/$type/g;

    return @count;
}

sub make_fasta_iterator {
    my ($fasta_file) = @_;

    open my $FASTA, '<', $fasta_file or croak "Can't read $fasta_file: $!";

    return sub {
        my $header   = <$FASTA>;
        my $sequence = <$FASTA>;

        defined $header   ? chomp $header   : return;
        defined $sequence ? chomp $sequence : return;

        ($header =~ s/>//) 
        or croak "Something wrong with fasta record: header $header does not seem to be a valid fasta header";

        return {header => $header, sequence => $sequence};
    };

}


__END__


=head1 NAME

 fasta_motifs.pl - Count frequencies of n motifs in m fasta files

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

 fasta_motifs.pl [OPTION]... [FILE]...

 -m, --motifs      any number of strings of motifs separated by spaces
 -f, --fasta       any number of fasta files
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
