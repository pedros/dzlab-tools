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


while (my $alignments = next_multiple_alignments( $left ) ) {
    #die Dumper $alignments->[0]->strand;
    print Dumper $alignments->[0]->strand;
}





    ##### No matches on either end #####

    ##### Unique matches on both ends #####

    ##### No match on one end, 1 match on other #####

    ##### No match on one end, multiple matches on other #####

    ##### One match on one end, multiple matches on other end #####

    ##### Multiple matches on both ends #####


=head2 choose_mates

 Takes an iterator of two Bio::DB::Bam::AlignWrapper objects as produced by make_mates_iterator(),
 and an options hash with the following keys:

 %opt = (
     random => boolean,
     insert-size => integer,
     variance    => float coefficient (0.0-1.0)
     orientation => -1|1,
 );

 Returns best pair if unique possible mapping,
 a random pair if multiple possible mappings and $opt{random},
 or undef if no possible mappings exist.

=cut

sub choose_mates {
    my ($alignment_iterator, %opts) = @_;

    my @alignments;

    while ( my $mates = $alignment_iterator->() ) {
        push @alignments, $mates
        if check_mates( @$mates, %opts );
    }
    
    return shift @alignments if @alignments == 1;
    return $alignments[int rand (@alignments - 1)] if @alignments > 1 and $opts{random};
    return unless @alignments;
}


=head2 check_mates

 Takes two Bio::DB::Sam::AlignWrapper objects,
 and an options hash with the following keys:

 %opt = (
     insert-size => integer,
     variance    => float coefficient (0.0-1.0)
     orientation => -1|1,
 );

 Returns true if the Bio::DB::Sam::AlignWrapper objects:

 1. Map to same sequence ID (chromosome)
 2. Map to opposite (-1) or same (1) orientation
 3. Are separated, from end of one object to start of another,
    by exactly $insert-size bp * ( +/- $variance )

=cut

sub check_mates {
    my ($left, $right, %opts) = @_;

    return if $left->get_tag_values('UNMAPPED') or $right->get_tag_values('UNMAPPED');

    return unless $left->seq_id eq $right->seq_id                        # condition 1, map to same sequence
    and           $left->strand *  $right->strand == $opts{orientation}; # condition 2, map according to $opts{orientation}

    # HERE BE DRAGONS, breaking encapsulation when accessing {sam}.
    # accessing $left or $right is the same after the check above.
    my $chr_len = $left->{sam}->length($left->seq_id);

    my $insert = $opts{orientation} == 1
                 ? $right->start - $left->end
                 : abs( $chr_len - $right->start + 1 ) - $left->end;

    my $variance = $opts{distance} * $opts{variance};

    return ($opts{distance} + $variance) >= $insert
    and    ($opts{distance} - $variance) <= $insert; # condition 3, map within acceptable range
}

=head2 make_mates_iterator

 Takes two references to arrays of Bio::DB::Bam::AlignWrapper objects.
 Returns an anonymous subroutine that will iterate over all combinations
 of one object from each array reference.

=cut

sub make_mates_iterator {
    my ($left, $right) = @_;

    my @combos;

    while (my $l = shift @$left) {
        for my $r (@$right) {
            push @combos, [$l, $r];
        }
    }

    return sub { shift @combos };
}


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
                my @reads        = @{$buffer{$sam}};
                @{$buffer{$sam}} = $read;
                return \@reads
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
