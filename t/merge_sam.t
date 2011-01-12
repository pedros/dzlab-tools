#!/usr/bin/env perl

use strict;  use warnings;
use Test::More qw/no_plan/;
use FindBin; use lib "$FindBin::Bin/..";

use Data::Dumper;

require_ok( 'merge_sam.pl' );

my ($left_sam, $right_sam) = map {
    Bio::DB::Sam->new( 
        -bam => $_,
        -expand_flags => 1,
        #-fasta => $ARGV{reference},
    )
} glob( "$FindBin::Bin/data/*bam" );

my ($left_sam_it, $right_sam_it) = map {
    $_->features( -iterator => 1 )
} ($left_sam, $right_sam);

my ($left_test_data, $right_test_data) = map {
    my $file = $_;
    [
        map { chomp; [split /\t/] }
        qx[ $FindBin::Bin/data/process_sam.pl $file ]
     ]
 } glob( "$FindBin::Bin/data/*sam" );

test_next_multiple_alignments( $left_sam_it, @$left_test_data );
test_next_multiple_alignments( $right_sam_it, @$right_test_data );

sub test_next_multiple_alignments {
    my ($sam_it, @test_data) = @_;

    while (my $alignments = merge_sam::next_multiple_alignments( $sam_it ) ) {

        my ($expected_seqid, $expected_count) = @{ shift @test_data };

        for my $alignment ( @$alignments ) {
            is( $alignment->query->seq_id, $expected_seqid, "Correct aggregate seqid" );
        }

        is( $expected_count, @$alignments, "Correct aggregate seqid count" );
    }
}

#die Dumper $right_test_data;


#merge_sam::choose_mates();
