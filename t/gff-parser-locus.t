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
use GFF::Parser::Locus;

my $p = GFF::Parser::Locus->new(file => \*DATA, locus => 'ID');

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
            attribute => "ID= mRNA0001.3 ;  c  =  123  ;  n  =  1234  ;  t  =  4  ",
            c => "123",
            n => 1234,
            t => "4",
            ID => "mRNA0001.3",
            'ID.prefix' => "mRNA0001",
            'ID.suffix' => "3",
        }, "gff line 1- trailing && embedded whitespace");
}

__DATA__
hello	.	exon	199	1233	1.2	-	2	ID= mRNA0001.3 ;  c  =  123  ;  n  =  1234  ;  t  =  4  
