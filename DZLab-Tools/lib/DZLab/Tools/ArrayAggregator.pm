package ArrayAggregator;

=head1 METHODS

=head2 new

Return new empty blessed []

=cut
sub new { bless [], shift }

=head2 step @values

Push stringified @values to $self

=cut
sub step {
    my ( $self, @values ) = @_;
  
    push @$self, "@values";
}


=head2 finalize

Stringify all elements in $self. This needs to happen because finalize is expected to return an aggregate, ie. a scalar value.
As such, list context would be coerced to the number of elements return, and an arrayref would be stringified.

=cut
sub finalize {
    my ($self) = @_;

    join "\t", @$self;
}


=head2 post_process $aggregate $by $cols

$aggregate is a row_hashref as returned by DBI. 
$by is the group-by col set in the query. 
$cols is the comma-separated stringified list of cols to aggregate on set in the query.

Returns [ranges ...], $by

=cut
sub post_process {
    my ($self, $aggregate, $by, $cols) = @_;

    my @ranges = split /\t/, $aggregate->{"aggregate($cols)"};
    @ranges = map { [split / /, $_] } @ranges;
    return \@ranges, $aggregate->{$by};
}

1;
