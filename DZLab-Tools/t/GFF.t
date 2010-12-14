#!/usr/bin/env perl

use FindBin;
use lib "$FindBin::Bin/../lib";
use DZLab::Tools::GFF;

use Test::Exception;
use Test::More;

BEGIN { use_ok( 'DZLab::Tools::GFF' ); }
require_ok( 'DZLab::Tools::GFF' );

my @gff_read_tests = qw/
                 blank blank_with_space
                 empty_comment empty_comment_with_prefixed_space empty_comment_with_suffixed_space
                 comment comment_with_prefixed_space comment_with_suffixed_space
                 empty_pragma empty_pragma_with_prefixed_space empty_pragma_with_suffixed_space
                 single_value_pragma multiple_value_pragma
                 gff_data_line
             /;

my %gff_read_common_struct = (
    source => q{.}, # TODO: substitute undef for q{.} once GFF.pm returns undefs,
    frame => q{.},
    score => q{.},
    end => 9000,
    feature => 'mRNA',
    seqname => 'ctg123',
    strand => '+',
    start => '1050'
);

my @gff_read_structs = (
    [], [], [], [], [], [], [], [], [], [], [],
    [qw/gff-version 3/], [qw/sequence-region ctg123 1 1497228/],
    {
        %gff_read_common_struct,
        attribute => 'ID=mRNA00001;Parent=gene00001;Name=EDEN.1',
        attributes => {
            ID => 'mRNA00001',
            Name => 'EDEN.1',
            Parent => 'gene00001'
        },
    },
    {
        %gff_read_common_struct,
        attribute => 'ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123',
        attributes => {
            ID => 'mRNA00001',
            Name => ['EDEN.1', 123],
            Parent => 'gene00001'
        },
    },
    {
        %gff_read_common_struct,
        attribute => 'mRNA00001',
        attributes => {Note => 'mRNA00001'},
    }
);

my @gff_make_iterator_bad_tests = qw/
                                        no_options
                                        no_sub_as_parser
                                        no_file_as_file
                                        no_glob_as_handle
                                        sub_as_parser_but_no_file_or_handle
                                    /;


my @gff_make_iterator_good_tests = qw/
                                        parser_and_file_ok
                                        parser_and_file_ok_bad_handle
                                        parser_and_handle_ok
                                        parser_and_handle_ok_bad_file
                                        parser_and_handle_and_file_ok
                                    /;

my @gff_make_iterator_bad_options = (
    {},
    {parser => 'not a function'},
    {file => 'not a file'},
    {handle => 'not a handle'},
    {parser => \&gff_read},
);

my @gff_make_iterator_good_options = (
    {parser => \&gff_read, file => $0},
    {parser => \&gff_read, file => $0, handle => 'not a handle'},
    {parser => \&gff_read, handle => \*DATA},
    {parser => \&gff_read, handle => \*DATA, file => 'not a file'},
    {parser => \&gff_read, handle => \*DATA, file => $0},
);


throws_ok(
    sub { gff_make_iterator( %{ shift @gff_make_iterator_bad_options } ) },
    qr/Need parser function reference and file name or handle to build iterator/,
    shift @gff_make_iterator_bad_tests
)
while @gff_make_iterator_bad_tests;

is(
    ref gff_make_iterator( %{ shift @gff_make_iterator_good_options } ),
    'CODE',
    shift @gff_make_iterator_good_tests
)
while @gff_make_iterator_good_tests;

my $it = gff_make_iterator(
    parser => \&gff_read,
    handle => \*DATA,
);

while (my $gff_struct = $it->()) {
    is_deeply(
        $gff_struct,
        shift @gff_read_structs, 
        shift @gff_read_tests
    )
}


done_testing();


__DATA__

 
#
 #
 # 
#helo
 #helo
 #helo asd 
##
 ##
 ## 
##gff-version 3
  ## sequence-region   ctg123 1 1497228
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123
ctg123	.	mRNA	1050	9000	.	+	.	mRNA00001
