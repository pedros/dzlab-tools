#!/usr/bin/env perl

use warnings;     use strict; use diagnostics;
use Data::Dumper; use Carp;
use Getopt::Long; use Pod::Usage;
use version;      our $VERSION = qv('0.0.1');
use FindBin;      use lib "$FindBin::Bin/DZLab-Tools/lib";

use DZLab::Tools::RunUtils;
use Bio::DB::Sam;

GetOptions(
    \%ARGV,
    'forward-reverse|fr', 'forward-forward|ff',
    'insert-size|i=i',    'reference|r=s',
    'output|o=s',         'error|e=s',
    _meta_options( \%ARGV ),
) and ($ARGV[0] and $ARGV[1] ) or pod2usage( -verbose => 1 );

my ( $INH, $OUTH, $ERRH ) = _prepare_io( \%ARGV, \@ARGV );

my ($left, $right) = map {
    Bio::DB::Sam->new( 
        -bam => $_,
        -expand_flags => 1,
        -fasta => $ARGV{reference},
    )->features( -iterator => 1 ) 
} @ARGV[0, 1];

while (my @alignments = next_multiple_alignments( $left ) ) {

}




    ##### No matches on either end #####

    ##### Unique matches on both ends #####

    ##### No match on one end, 1 match on other #####

    ##### No match on one end, multiple matches on other #####

    ##### One match on one end, multiple matches on other end #####

    ##### Multiple matches on both ends #####



=head2 next_multiple_alignments

 Takes a Bio::DB::Sam features iterator.
 Returns undef if iterator exhausted,
 or a list of Bio::DB::Bam::AlignWrapper objects
 representing different possible Sam alignments
 for the same read (as defined by its group/machine id)

=cut
sub next_multiple_alignments {
    my ($sam) = @_;

    # `buffer` is of the form:
    # $buffer{$stringified_sam_iterator}[Bio::DB::Bam::AlignWrapper ...]
    use feature 'state';
    state %buffer;

    while ( my $read = $sam->next_seq ) {

        # not seen this iterator yet, or buffer empty
        # save read; note, buffer is stateful, so 
        # read will be possibly be available next call
        if (not exists $buffer{$sam} or not @{$buffer{$sam}}) {
            push @{$buffer{$sam}}, $read;
        }

        # buffer not empty: compared previous and current read ids
        # if they match, add current for later
        # if not, flush buffer, save current for later, return
        else {
            my $previous_read_id = $buffer{$sam}[-1]->query->seq_id;
            my $current_read_id  = $read->query->seq_id;

            if ($previous_read_id eq $current_read_id) {
                push @{$buffer{$sam}}, $read;
            }
            else {
                my @reads         = @{$buffer{$sam}};
                @{$buffer{$sam}} = $read;
                return @reads
            }
        }
    }
}




__END__

=head1 NAME

 APerlyName.pl - Short description

=head1 SYNOPSIS

 APerlyName.pl [OPTION]... [[-i] FILE]...

=head1 DESCRIPTION

 Long description

=head1 OPTIONS

 -i, --input       <string>  (STDIN)   input filename
 -o, --output      <string>  (STDOUT)  output filename
 -e, --error       <string>  (STDERR)  output error filename
     --verbose     [integer] (0)       print increasingly verbose error messages
     --quiet                           print no diagnostic or warning messages
     --debug                           run in debug mode
     --version                         print current version
     --license                         print author's contact and copyright information
     --help                            print this information
     --manual                          print the plain old documentation page

=head1 VERSION

 0.0.1

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
