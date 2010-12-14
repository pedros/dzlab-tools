package DZLab::Tools::RangeUtils;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use List::Util qw/max min/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(
overlap_ratio range_overlap uniq_ranges
range_before range_after
);

=head1 uniq_ranges [[$start_a, $end_a], [$start_b, $end_b], ... ]

Gets rid of duplicates from range list.

=cut 
sub uniq_ranges {
    my ($ranges) = @_;

    my %seen;
    my @uniq;

    for (@$ranges) {
        push @uniq, $_ unless $seen{"$_->[0];$_->[1]"}++;
    }

    return wantarray ? @uniq : [@uniq];
}

=head1 range_before [$start_a, $end_a], [$start_b, $end_b]

Return true if first range completley before second range

=cut 
sub range_before{
    my ($range_a, $range_b) = @_;

    return $range_a->[1] < $range_b->[0];
}

=head1 range_after [$start_a, $end_a], [$start_b, $end_b]

Return true if first range completley after second range

=cut 
sub range_after{
    my ($range_a, $range_b) = @_;

    return $range_a->[0] > $range_b->[1];
}

=head1 range_overlap [$start_a, $end_a], [$start_b, $end_b]

Return true if range_overlap shares even a single position.

=cut
sub range_overlap{
    my ($range_a, $range_b) = @_;

    return ! (
        range_before($range_a, $range_b) ||
        range_after($range_a, $range_b));
}

=head1 overlap_ratio [$start_a, $end_a], [$start_b, $end_b]

    a: |--------------|
    b:    |----------------------------|
          |-----------| <- overlap region

return ratio of overlap region length / range_b length, or 0 if no overlap

=cut 
sub overlap_ratio {
    
    my ($range_a, $range_b) = @_;
    return 0 if ! range_overlap($range_a,$range_b);

    my $bin_low  = max( $range_a->[0], $range_b->[0] );
    my $bin_high = min( $range_a->[1],   $range_b->[1]   );

    my $overlap  = abs($bin_high - $bin_low) + 1;

    return $overlap / ( $range_b->[1] - $range_b->[0] + 1 );
}

1;

