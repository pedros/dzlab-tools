#!/usr/bin/env perl
package GFF::Parser;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;

has filename_or_handle => (
    is => 'ro',
    required => 1,
    init_arg => 'file',
);

has locus_tag => (
    is => 'ro',
    init_arg => 'locus',
);

# privates

has filehandle => (
    is => 'rw',
    init_arg => undef,
);

sub BUILD{
    my ($self) = @_;
    if (ref $self->filename_or_handle eq 'GLOB'){
        $self->filehandle = $self->filename_or_handle;
    }
    elsif (!ref $self->filename_or_handle && -f $self->filename_or_handle ){
        open my $fh, '<', $self->filename_or_handle
            or croak "cannot open $self->filename_or_handle";
        $self->filehandle($fh);
    }
    else {
        croak "file argument to GFF::Parser needs to be file handle or file name";
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
        my $gff = _parse_gff_hashref($line,$self->locus_tag);
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
    while (defined (my $line = scalar readline $self->filehandle)){
        return _parse_gff_hashref($line,$self->locus_tag);
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
            if (exists $gff->{$column} && defined $gff->{$column}){
                push @{$index{$gff->{$column}}}, $gff;
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

=head2 parse_gff_hashref $line

read a line, and return a gff hashref with all attributes. attributes are returned the same level,
not in a sub hash.  maps '.' dot columns and non-existenct attributes to undef
returns 0 on non-parseable lines or comments. *Returns pragmas as strings*.
Check that ref $result eq 'HASH'

If a locus_tag is given as a second argument, and locus is of the form "SOME.THING",
locus_tag.prefix is set to SOME and locus_tag.suffix is set to THING.  If locus does 
not have a dot, locus_tag.prefix and suffix are undef.  

=cut

sub _parse_gff_hashref{
    my ($line, $locus_tag) = @_;
     
    return 0 unless $line;
    $line =~ tr/\n\r//d;

    if ($line =~ m/^\s*#(#)?/o){
        return defined $1 ? $' : 0; # $' = regex postmatch
    }

    my %accum;
    @accum{qw/seqname source feature start end score strand frame attribute/}
    = map { $_ eq q{.} ? undef : $_ } split /\t/, $line;

    (carp "unparseable GFF line" && return 0) unless keys %accum == 9;

    if (defined $accum{seqname}){
        $accum{seqname} = lc $accum{seqname};
    }

    if (defined($accum{attribute}) ){
        for (split /;/, $accum{attribute}){
            my ($key, $val) = split /=/, $_;
            $key =~ s/^\s+//;
            $key =~ s/\s+$//;

            if (defined $val){
                $val =~ s/^\s+//;
                $val =~ s/\s+$//;
                $accum{$key} = $val;
            }
            else {
                $accum{Note} = $key;
            }
        }
    }

    if ($locus_tag){
        if (exists($accum{$locus_tag}) && 
            $accum{$locus_tag} =~ /([^\.]+)\.([^\.]+)/){
            $accum{$locus_tag . '.prefix'} = $1;
            $accum{$locus_tag . '.suffix'} = $2;
        } else {
            $accum{$locus_tag . '.prefix'} = undef;
            $accum{$locus_tag . '.suffix'} = undef;
        }
    }

    return \%accum;
}

sub _is_gff{
    my ($gff) = @_;
    return ref $gff && ref $gff eq 'HASH';
}

no Moose;
__PACKAGE__->meta->make_immutable;

1;

