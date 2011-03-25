#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;

use Test::More qw(no_plan);

use FindBin;
use lib "$FindBin::Bin/../lib";

use GFF;
use GFF::Parser;

my $p = GFF::Parser->new(file => \*DATA);

### first 8 lines should be blank

for (1 .. 8){ 
    my $gff = $p->next_no_skip();
    is($gff,0, "comment/blank should be zero");
}

### pragmas

{
    my $gff = $p->next_no_skip();
    like($gff,qr/pragma 1/, "pragma test 1");
}

{
    my $gff = $p->next_no_skip();
    like($gff,qr/pragma 2/, "pragma test 2");
}

{
    my $gff = $p->next_no_skip();
    like($gff,qr/\s+/, "pragma test blank");
}

### no attirbutes. 

{
    my $gff = $p->next_no_skip();
    ok(!$gff->equals(
            sequence => 'ctg123-meow', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            attribute_string => "ID=mRNA00001;Parent=gene00001;Name=EDEN.1"
        ),
        "gff line 1- should NOT be equal");
}

{
    my $gff = $p->next_no_skip();
    ok($gff->equals(
            sequence => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            attribute_string => "ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123",
        )
        , "gff line 2- no attr");
}

{
    my $gff = $p->next_no_skip();
    ok($gff->equals(
            sequence => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            attribute_string => "mRNA00001",
        ), "gff line 3- no attr");
}

### With attributes. don't pass attribute_string b/c don't care if we have parsed attrs

{
    my $gff = $p->next_no_skip();
    ok($gff->equals(
            sequence => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            ID => "mRNA00001",
            Parent => "gene00001",
            Name => "EDEN.1",
        ), "gff with  attributes line 1");
}

{
    my $gff = $p->next_no_skip();
    ok($gff->equals(
            sequence => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            ID => "mRNA00001",
            Parent => "gene00001",
            Name => "EDEN.1,123",
        ), "gff with  attributes line 2");
}

{
    my $gff = $p->next_no_skip();
    ok($gff->equals(
            sequence => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            Note => "mRNA00001",
        ), "gff with  attributes line 3");
}

{
    my $gff = $p->next_no_skip();
    ok(!$gff->equals(
            sequence => 'hello', 
            source => undef, 
            feature => 'exon', 
            start => 199, 
            end => 1233, 
            frame => 2,
            strand => '-',
            score => 1.2,
            c => 123,
            n => 666,
            t => 4,
        ), "gff with  attributes line 4- should NOT be equal");
}

{
    my $gff = $p->next_no_skip();
    ok($gff->equals(
            sequence => 'hello', 
            source => undef, 
            feature => 'exon', 
            start => 199, 
            end => 1233, 
            frame => 2,
            strand => '-',
            score => 1.2,
            c => 123,
            n => 1234,
            t => 4,
        ), "gff with  attributes line 5- trailing whitespace");
}

{
    my $gff = $p->next_no_skip();
    ok($gff->equals(
            sequence => 'hello', 
            source => undef, 
            feature => 'exon', 
            start => 199, 
            end => 1233, 
            frame => 2,
            strand => '-',
            score => 1.2,
            c => "123",
            n => 1234,
            t => "4",
        ), "gff with  attributes line 5- trailing && embedded whitespace");
}

{
    my $gff = $p->next_no_skip();
    is_deeply([$gff->parse_locus('ID')],['mRNA0001','3'], "parsing locus");
}



__DATA__

 
#
 #
 # 
#helo
 #helo
 #helo asd 
## pragma 1
 ## pragma 2
 ## 
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123
ctg123	.	mRNA	1050	9000	.	+	.	mRNA00001
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123
ctg123	.	mRNA	1050	9000	.	+	.	mRNA00001
hello	.	exon	199	1233	1.2	-	2	c=123;n=1234;t=4
hello	.	exon	199	1233	1.2	-	2	c=123;n=1234;t=4  
hello	.	exon	199	1233	1.2	-	2	  c  =  123  ;  n  =  1234  ;  t  =  4  
hello	.	exon	199	1233	1.2	-	2	ID= mRNA0001.3 ;  c  =  123  ;  n  =  1234  ;  t  =  4  
