#!/usr/bin/env perl
package GFF::Parser::Attributes;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;

extends 'GFF::Parser';

=head2 parse_gff_hashref $line

read a line, and return a gff hashref with all attributes. attributes are returned the same level,
not in a sub hash.  maps '.' dot columns and non-existenct attributes to undef
returns 0 on non-parseable lines or comments. *Returns pragmas as strings*.
Check that ref $result eq 'HASH'

If a locus_tag is given as a second argument, and locus is of the form "SOME.THING",
locus_tag.prefix is set to SOME and locus_tag.suffix is set to THING.  If locus does 
not have a dot, locus_tag.prefix and suffix are undef.  

=cut

around '_parse_gff' => sub{
    my $orig = shift;
    my $self = shift;
    my $line = shift;

    my $gff = $self->$orig($line);
    
    if (ref $gff eq 'HASH'){
        if (defined($gff->{attribute}) ){
            for (split /;/, $gff->{attribute}){
                my ($key, $val) = split /=/, $_;
                $key =~ s/^\s+//;
                $key =~ s/\s+$//;

                if (defined $val){
                    $val =~ s/^\s+//;
                    $val =~ s/\s+$//;
                    $gff->{$key} = $val;
                }
                else {
                    $gff->{Note} = $key;
                }
            }
        }
    }

    return $gff;
};

no Moose;
__PACKAGE__->meta->make_immutable;

1;

