#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';

#use Test::Simple qw(no_plan);
use Test::More tests => 3;
BEGIN { use_ok( 'DZLab::Tools::GFF' ); }
require_ok( 'DZLab::Tools::GFF' );


my %gff_read_common_struct = (
    source => q{.}, # TODO: substitute undef for q{.} once GFF.pm returns undefs,
    frame => q{.},
    score => q{.},
    strand => '+',
);

my $expected = 
{ ctg123 => 
    [
    {
        %gff_read_common_struct,
        feature => 'mRNA',
        seqname => 'ctg123',
        end => 9000,
        start => 1050,
        range => [1050,9000],
        attribute => 'ID=mRNA00001;Parent=gene00001;Name=EDEN.1',
        attributes => {
            ID => 'mRNA00001',
            Name => 'EDEN.1',
            Parent => 'gene00001'
        },
    },
    {
        %gff_read_common_struct,
        feature => 'mRNA',
        seqname => 'ctg123',
        end => 9000,
        start => 1051,
        range => [1051,9000],
        attribute => 'ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123',
        attributes => {
            ID => 'mRNA00001',
            Name => ['EDEN.1', 123],
            Parent => 'gene00001'
        },
    },
    {
        %gff_read_common_struct,
        feature => 'mRNA',
        seqname => 'ctg123',
        end => 9000,
        start => 1052,
        range => [1052,9000],
        attribute => 'mRNA00001',
        attributes => {Note => 'mRNA00001'},
    },
    ],
    foobar => 
    [
    {
        %gff_read_common_struct,
        end => 9000,
        feature => 'codon',
        seqname => 'foobar',
        start => 1123,
        range => [1123,9000],
        attribute => 'ID=mRNA00002',
        attributes => {ID => 'mRNA00002'},
    },
    {
        %gff_read_common_struct,
        end => 9000,
        feature => 'codon',
        seqname => 'foobar',
        start => 1251,
        range => [1251,9000],
        attribute => 'ID=mRNA00002',
        attributes => {ID => 'mRNA00002'},
    },
    {
        %gff_read_common_struct,
        end => 9000,
        feature => 'codon',
        seqname => 'foobar',
        start => 1459,
        range => [1459,9000],
        attribute => 'mRNA00002',
        attributes => {Note => 'mRNA00002'},
    },
    ]
};


my $read = gff_slurp_by_seq {handle => \*DATA, debug => 1,sort => 1};

is_deeply($read,$expected);


__DATA__
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123	.	mRNA	1051	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123
ctg123	.	mRNA	1052	9000	.	+	.	mRNA00001
foobar	.	codon	1251	9000	.	+	.	ID=mRNA00002
foobar	.	codon	1123	9000	.	+	.	ID=mRNA00002
foobar	.	codon	1459	9000	.	+	.	mRNA00002
