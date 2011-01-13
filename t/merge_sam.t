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

my $left_alignments  = test_next_multiple_alignments( $left_sam_it, @$left_test_data );
my $right_alignments = test_next_multiple_alignments( $right_sam_it, @$right_test_data );
my $mates_iterators  = test_make_mates_iterator( $left_alignments, $right_alignments );


test_check_mates( $mates_iterators );


sub test_check_mates {
    my ($mates_iterators) = @_;

    my %opts = (
        insert   => 0,
        variance => 0,
        orientation => 1,
    );

    while (my $mates_iterator = shift @$mates_iterators ) {

        while (my $mates = $mates_iterator->() ) {

            if (merge_sam::check_mates( @$mates, %opts )) {
                is( $mates->[0]->strand, $mates->[1]->strand, "Accepted mates on same orientation" );
                is( $mates->[0]->end, $mates->[1]->start, "Accepted mates are adjacent" );
            }
        }
    }
}



sub test_make_mates_iterator {
    my ($left_alignments, $right_alignments) = @_;

    my @iterators;

    while( my $left_alignment = shift @$left_alignments and my $right_alignment = shift @$right_alignments ) {

        my $expected_number_of_alignments = @$left_alignment * @$right_alignment;
        my $actual_number_of_alignments   = 0;

        my $mates_it = merge_sam::make_mates_iterator( [@$left_alignment], [@$right_alignment] );

        while( my $mates = $mates_it->() ) {
            is( $mates->[0]->query->seq_id, $mates->[1]->query->seq_id, "Mates have same machine ID" );
            $actual_number_of_alignments++;
        }

        is( $expected_number_of_alignments, $actual_number_of_alignments, "Correct number of mate-mate combinations" );

        push @iterators, merge_sam::make_mates_iterator( $left_alignment, $right_alignment );
    }
    return \@iterators;
}


sub test_next_multiple_alignments {
    my ($sam_it, @test_data) = @_;

    my @alignments;

    while (my $alignments = merge_sam::next_multiple_alignments( $sam_it ) ) {

        my ($expected_seqid, $expected_count) = @{ shift @test_data };

        for my $alignment ( @$alignments ) {
            is( $alignment->query->seq_id, $expected_seqid, "Correct aggregate seqid" );
        }

        is( $expected_count, @$alignments, "Correct aggregate seqid count" );

        push @alignments, $alignments;
    }
    return \@alignments;
}

#die Dumper $right_test_data;


#merge_sam::choose_mates();

