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
has sequence => (is => 'ro',isa => 'Str');
has window_args => (is => 'ro',default => sub { [] });

# return a Window (or descendant) for next one in line
sub next{ return; }
sub create_window{
    my ($self,$start,$end,$window_id)=@_;
    return $self->class->new(
        sequence => $self->sequence, 
        start => $start,
        end => $end,
        window_id => $window_id,
        @{$self->window_args});
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

#==================================================================
# WindowSource::Fixed

=head2 WindowSource::Fixed

    my $ws = WindowSource::Fixed->new(
        class => 'Window::SumScore', 
        length => 1000, 
        step => 123, 
        sequence => 'chr1',
    );
    while (my $window = $ws->next()){
        say $window->to_gff;
    }

=cut

package WindowSource::Fixed;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use autodie;
extends 'WindowSource';

has current => (
    traits  => ['Counter'],
    is => 'rw', 
    isa => 'Int', 
    default => 1,
    handles => {
        inc_current   => 'inc',
    },
    init_arg => undef,
);
has length  => (is => 'ro', isa => 'Int');
has step    => (is => 'ro', isa => 'Int');

override 'next' => sub{
    my ($self) = @_;
    my $start = $self->current;

    return if $self->current > $self->length;

    my $end = $self->current + $self->step - 1;

    if ($end > $self->length){
        $end = $self->length;
    }
    $self->inc_current($self->step);

    return $self->create_window($start,$end,"[$start,$end]");
};

no Moose;
__PACKAGE__->meta->make_immutable;

1;

#==================================================================
# WindowSource::Annotation

=head2 WindowSource::Annotation

    my $wsa = WindowSource::Annotation->new(
        class => 'Window::SumScore',
        file => 'tmp/TAIR8_genes-chrc.gff.sorted',
        locus => 'ID',
        sequence => 'chrc'
    );
    while (my $window = $wsa->next()){
        say $window->to_gff;
    }

=cut

package WindowSource::Annotation;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use autodie;
extends 'WindowSource';

has file => (is => 'ro',required => 1);
has parser => (is => 'rw', isa => 'GFF::Parser');
has locus => (is => 'ro', isa => 'Str');

override 'next' => sub{
    my ($self) = @_;
    my $gff = $self->parser->next();
    if (ref $gff eq 'GFF'){
        return $self->create_window($gff->start,$gff->end,$gff->get_column($self->locus));
    }
};
sub BUILD{
    my $self=shift;
    $self->parser(GFF::Parser->new(file => $self->file));
}

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
use Fasta;

use Window;

unless (caller){
    my $ws = WindowSource::Fixed->new(
        class => 'Window::SumScore', 
        length => 1000, 
        step => 123, 
        sequence => 'chr1',
    );
    while (my $window = $ws->next()){
        say $window->to_gff;
    }
    my $wsa = WindowSource::Annotation->new(
        class => 'Window::SumScore',
        file => 'tmp/TAIR8_genes-chrc.gff.sorted',
        locus => 'ID',
        sequence => 'chrc'
    );
    while (my $window = $wsa->next()){
        say $window->to_gff;
    }

        
    
}


