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
    print Dumper $self;
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

sub next{
    my ($self) = @_;
    my $line = scalar readline $self->filehandle;
    if (defined ($line)){
        return _parse_gff_hashref($line,$self->locus_tag);
    }
    else {
        return;
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
            if (defined $val){
                $accum{$key} = $val;
            }
            else {
                $accum{Note} = $key;
            }
        }
    }

    if ($locus_tag){
        if ($accum{$locus_tag} =~ /([^\.]+)\.([^\.]+)/){
            $accum{$locus_tag . '.prefix'} = $1;
            $accum{$locus_tag . '.suffix'} = $2;
        } else {
            $accum{$locus_tag . '.prefix'} = undef;
            $accum{$locus_tag . '.suffix'} = undef;
        }
    }

    return \%accum;
}


no Moose;
__PACKAGE__->meta->make_immutable;

1;

