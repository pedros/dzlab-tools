#!/usr/bin/env perl
package GFF::Parser;
use Moose;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use GFF;

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
        my $gff = $self->_parse_gff($line);
        if (_is_gff($gff)){
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
        return $self->_parse_gff($line);
    }
    return;
}

=head2 $p->slurp()

return an arrayref of all gff records

=cut

sub slurp{
    my ($self) = @_;
    my @accum;
    while (my $gff = $self->next()){
        # if we're skipping, no need to check
        if (_is_gff($gff)){ 
            push @accum, $gff;
        }
    }
    return \@accum;
}

=head2 $p->slurp_index('colname')

return a hashref of column val to gff record

=cut

sub slurp_index{
    my ($self, $column) = @_;
    my %index;
    my $counter = 0;
    my $badrecords = 0;
    while (my $gff = $self->next()){
        if (_is_gff($gff)){ 
            if (defined $gff->get_column($column)){
                push @{$index{$gff->get_column($column)}}, $gff;
            } else {
                $badrecords++;
            }
            $counter++;
        }
    }
    carp "warning: $badrecords out of $counter records didn't have a column/attribute $column in " .  $self->filename_or_handle
    if $badrecords;

    return \%index;
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

=head2 $parser->_parse_gff($line)

read a line, and return a gff hashref.  maps '.' dot columns and non-existenct attributes to undef. returns 0 on
non-parseable lines or comments. *Returns pragmas as strings*.  Check that ref $result eq 'HASH'

it's a member method so it can be inherited

=cut

sub _parse_gff{
    my ($self,$line) = @_;
     
    return 0 unless $line;
    $line =~ tr/\n\r//d;

    if ($line =~ m/^\s*#(#)?/o){
        return defined $1 ? $' : 0; # $' = regex postmatch
    }

    my @split = split /\t/, $line;
    (carp "unparseable GFF line: $line" && return 0) unless @split == 9;

    my %accum;
    @accum{qw/sequence source feature start end score strand frame attribute_string/}
    = map { $_ eq q{.} ? undef : $_ } @split;

    if (defined $accum{sequence}){
        $accum{sequence} = lc $accum{sequence};
    }

    return GFF->new(%accum);
}

sub _is_gff{
    my ($gff) = @_;
    return ref $gff eq 'GFF';
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

