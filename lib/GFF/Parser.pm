#!/usr/bin/env perl
package GFF::Parser;
use Moose;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use GFF;
use GFF::Util;
use autodie;

has filename_or_handle => (
    is => 'ro',
    required => 1,
    init_arg => 'file',
);

# privates

has filehandle => (
    is => 'rw',
    init_arg => undef,
);

sub BUILD{
    my ($self) = @_;
    if (ref $self->filename_or_handle eq 'GLOB'){
        $self->filehandle($self->filename_or_handle);
    }
    elsif (!ref $self->filename_or_handle && -f $self->filename_or_handle ){
        open my $fh, '<', $self->filename_or_handle
            or croak "cannot open $self->filename_or_handle";
        $self->filehandle($fh);
    } elsif (! -f $self->filename_or_handle){
        croak $self->filename_or_handle . " doesn't exist?";
    }
    else {
        croak "file argument to GFF::Parser needs to be file handle or file name" 
        . Dumper $self;
    }
}

sub DEMOLISH{
    my ($self) = @_;
    if (!ref $self->filename_or_handle){
        close $self->filehandle
            or croak "cannot close $self->filename_or_handle";
    }
}

=head2 $p->next()

return the next gff record. returns undef when eof.

=cut 

sub next{
    my ($self) = @_;
    while (defined (my $line = scalar readline $self->filehandle)){
        my $gff = parse_gff($line);
        if (is_gff($gff)){
            return $gff;
        }
    }
    return;
}

=head2 $p->next_no_skip()

return the next gff record, 0 on a comment, or string on a pragma statement.
returns undef when eof.

=cut 

sub next_no_skip{
    my ($self) = @_;
    while (defined (my $line = readline $self->filehandle)){
        return parse_gff($line);
    }
    return;
}

=head2 $p->do(sub { my $gff = shift; ... })

=cut

sub do{
    my ($self,$code) = @_;
    croak "do needs a sub!" if (ref $code ne 'CODE');
    while (my $gff = $self->next()){
        $code->($gff);
    }
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

