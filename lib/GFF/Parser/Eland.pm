#!/usr/bin/env perl
package GFF::Parser::Eland;
use Moose::Role;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';

has locus_tag => (
    is => 'ro',
    init_arg => 'locus',
);

around '_parse_gff' => sub{
    my $orig = shift;
    my $self = shift;
    my $line = shift;
    my $gff = $orig->$self($line);

    if ( ref $gff eq 'HASH'){
        $gff->{feature} =~ m{^([^/]+)/(\d+):[natgcNATGC]+$};
        @{$gff}{'readid','pairid','readseq'} = ($1,$2,$3);

        #what is the source column?
    }
    return $gff;
};

1;


# chr5    R/R     GA2_0029_FC:6:1:1410:950#CCCCGA/1:NTAATAAATGATATTGTGGATTCTGATGGTAAGGTTGTATAATTG 8847396 8847440 1       -       2       target=CCCCAACAAATGATATTGTGGACCCTGATGGCAAGGTTGCATAAATGCG
