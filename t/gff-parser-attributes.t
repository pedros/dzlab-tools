#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Smart::Comments;

use Test::More qw(no_plan);

use FindBin;
use lib "$FindBin::Bin/../lib";

use GFF;
use GFF::Parser::Attributes;

my $p = GFF::Parser::Attributes->new(file => \*DATA);

{
    my $gff = $p->next_no_skip();
    is_deeply($gff,{
            seqname => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            attribute => "ID=mRNA00001;Parent=gene00001;Name=EDEN.1",
            ID => "mRNA00001",
            Parent => "gene00001",
            Name => "EDEN.1",
        }, "gff line 1");
}

{
    my $gff = $p->next_no_skip();
    is_deeply($gff,{
            seqname => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            attribute => "ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123",
            ID => "mRNA00001",
            Parent => "gene00001",
            Name => "EDEN.1,123",
        }, "gff line 2");
}

{
    my $gff = $p->next_no_skip();
    is_deeply($gff,{
            seqname => 'ctg123', 
            source => undef, 
            feature => 'mRNA', 
            start => 1050, 
            end => 9000, 
            frame => undef,
            strand => '+',
            score => undef,
            attribute => "mRNA00001",
            Note => "mRNA00001",
        }, "gff line 3");
}

{
    my $gff = $p->next_no_skip();
    is_deeply($gff,{
            seqname => 'hello', 
            source => undef, 
            feature => 'exon', 
            start => 199, 
            end => 1233, 
            frame => 2,
            strand => '-',
            score => 1.2,
            attribute => "c=123;n=1234;t=4",
            c => 123,
            n => 1234,
            t => 4,
        }, "gff line 4");
}

{
    my $gff = $p->next_no_skip();
    is_deeply($gff,{
            seqname => 'hello', 
            source => undef, 
            feature => 'exon', 
            start => 199, 
            end => 1233, 
            frame => 2,
            strand => '-',
            score => 1.2,
            attribute => "c=123;n=1234;t=4  ",
            c => 123,
            n => 1234,
            t => 4,
        }, "gff line 5- trailing whitespace");
}

{
    my $gff = $p->next_no_skip();
    is_deeply($gff,{
            seqname => 'hello', 
            source => undef, 
            feature => 'exon', 
            start => 199, 
            end => 1233, 
            frame => 2,
            strand => '-',
            score => 1.2,
            attribute => "  c  =  123  ;  n  =  1234  ;  t  =  4  ",
            c => "123",
            n => 1234,
            t => "4",
        }, "gff line 5- trailing && embedded whitespace");
}

__DATA__
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1,123
ctg123	.	mRNA	1050	9000	.	+	.	mRNA00001
hello	.	exon	199	1233	1.2	-	2	c=123;n=1234;t=4
hello	.	exon	199	1233	1.2	-	2	c=123;n=1234;t=4  
hello	.	exon	199	1233	1.2	-	2	  c  =  123  ;  n  =  1234  ;  t  =  4  
