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
use GFF::Parser;

my $p = GFF::Parser->new(file => \*DATA);

### first 8 lines should be blank

for (1 .. 8){ 
    my $gff = $p->next_no_skip();
    is($gff,0, "comment/blank should be zero");
}

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
        }, "gff line 3");
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
