#!/usr/bin/env perl
# support file for window_gff_new.pl

#==================================================================
# WindowSource

package WindowSource;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use autodie;

has class => (is => 'ro');
has seqname => (is => 'ro',isa => 'Str');

sub BUILD{
}

# return a Window (or descendant) for next one in line
sub next{
    my ($self) = @_;
    return;
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

#==================================================================
# WindowSource::Fixed

package WindowSource::Fixed;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use autodie;
extends 'WindowSource';

has 'current' => (is => 'rw', isa => 'Int', default => 1);
has 'length'  => (is => 'ro', isa => 'Int');
has 'step'    => (is => 'ro', isa => 'Int');

override 'next' => sub{
    my ($self) = @_;
    my $start = $self->current;

    return if $self->current > $self->length;

    my $end = $self->current + $self->step - 1;

    if ($end > $self->length){
        $end = $self->length;
    }
    $self->current($self->current+$self->step);

    return $self->class->new(seqname => $self->seqname, start => $start,end => $end);
};

no Moose;
__PACKAGE__->meta->make_immutable;

1;

#==================================================================
# WindowSource::Annotation

package WindowSource::Annotation;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use autodie;
extends 'WindowSource';

has 'file' => (is => 'ro');

override 'next' => sub{
    my ($self) = @_;
    #return $self->WindowClass->new(start => $start,end => $end);
    return;
};

no Moose;
__PACKAGE__->meta->make_immutable;

1;

#==================================================================
# main

package main;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

use FindBin;
use lib "$FindBin::Bin";

use Window;

unless (caller){
    my $ws = WindowSource::Fixed->new(class => 'Window::SumScore', length => 1000, step => 123, seqname => 'chr1');
    while (my $window = $ws->next()){
        say Dumper $window;
    }
}


