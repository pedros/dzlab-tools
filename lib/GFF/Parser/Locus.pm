#!/usr/bin/env perl
package GFF::Parser::Locus;
use Moose;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;

extends 'GFF::Parser::Attributes';

has locus_tag => (
    is => 'ro',
    init_arg => 'locus',
);

# extend _parse_gff to split the $locus_tag field into $locus_tag.prefix and $locus_tag.suffix, based on a separator
# (currently hard-coded to be a dot).
around '_parse_gff' => sub{
    my $orig = shift;
    my $self = shift;
    my $line = shift;
    my $gff = $self->$orig($line);
    my $locus_tag = $self->locus_tag;

    if ( ref $gff eq 'HASH'){
        if ($self->locus_tag()){
            if (exists($gff->{$locus_tag}) && 
                $gff->{$locus_tag} =~ /([^\.]+)\.([^\.]+)/){
                $gff->{$locus_tag . '.prefix'} = $1;
                $gff->{$locus_tag . '.suffix'} = $2;
            } else {
                $gff->{$locus_tag . '.prefix'} = undef;
                $gff->{$locus_tag . '.suffix'} = undef;
            }
        }
    }
    return $gff;
};

1;
