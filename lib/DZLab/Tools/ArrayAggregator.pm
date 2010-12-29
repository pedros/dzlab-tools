package ArrayAggregator;

sub new { bless [], shift }

sub step {
    my ( $self, @values ) = @_;
  
    push @$self, "@values";
}

sub finalize {
    my ($self) = @_;

    join "\t", @$self;
}

sub post_process {
    my ($self, $aggregate, $cols) = @_;

    my @ranges = split /\t/, $aggregate->{"aggregate($cols)"};
    @ranges = map { [split / /, $_] } @ranges;
    return [\@ranges, $aggregate->{$by}];
}

1;
