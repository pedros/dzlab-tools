#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use FindBin;
use lib "$FindBin::Bin/../lib";

use DZLab::Tools::GFFStore;

BEGIN { use_ok( 'DZLab::Tools::ArrayAggregator' ); }
require_ok( 'DZLab::Tools::ArrayAggregator' );

my $gff_store = DZLab::Tools::GFFStore->new({
    attributes => {
        Parent => 'text'
    },
    indices => [[qw/feature start end Parent/]],
});

$gff_store->slurp({ handle => \*DATA });

my $s = $gff_store->select_aggregate( 'exon', 'Parent', 'start', 'end' );

my %counts = (
    'AT1G01010.1' => 6,
    'AT1G01020.1' => 9,
    'AT1G01020.2' => 8,
    'AT1G01030.1' => 2,
    'AT1G01040.1' => 16,
);

my @ids = sort keys %counts;

while (my $locus = $s->()) {
    my ($ranges, $id) = @$locus;

    ok( exists $counts{$id},               'Exists aggregate ID'          );
    is( $id,                 shift @ids,   'Correct aggregate ID'         );
    is( scalar @$ranges,     $counts{$id}, 'Correct aggregated instances' );
};

done_testing();

__DATA__
Chr1	TAIR8	mRNA	3631	5899	.	+	.	ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Index=1
Chr1	TAIR8	protein	3760	5630	.	+	.	ID=AT1G01010.1-Protein;Name=AT1G01010.1;Derives_from=AT1G01010.1
Chr1	TAIR8	exon	3631	3913	.	+	.	Parent=AT1G01010.1
Chr1	TAIR8	five_prime_UTR	3631	3759	.	+	.	Parent=AT1G01010.1
Chr1	TAIR8	CDS	3760	3913	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR8	exon	3996	4276	.	+	.	Parent=AT1G01010.1
Chr1	TAIR8	CDS	3996	4276	.	+	2	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR8	exon	4486	4605	.	+	.	Parent=AT1G01010.1
Chr1	TAIR8	CDS	4486	4605	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR8	exon	4706	5095	.	+	.	Parent=AT1G01010.1
Chr1	TAIR8	CDS	4706	5095	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR8	exon	5174	5326	.	+	.	Parent=AT1G01010.1
Chr1	TAIR8	CDS	5174	5326	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR8	exon	5439	5899	.	+	.	Parent=AT1G01010.1
Chr1	TAIR8	CDS	5439	5630	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR8	three_prime_UTR	5631	5899	.	+	.	Parent=AT1G01010.1
Chr1	TAIR8	mRNA	6790	8737	.	-	.	ID=AT1G01020.1;Parent=AT1G01020;Name=AT1G01020.1;Index=1
Chr1	TAIR8	protein	6915	8666	.	-	.	ID=AT1G01020.1-Protein;Name=AT1G01020.1;Derives_from=AT1G01020.1
Chr1	TAIR8	five_prime_UTR	8667	8737	.	-	.	Parent=AT1G01020.1
Chr1	TAIR8	CDS	8571	8666	.	-	0	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR8	exon	8571	8737	.	-	.	Parent=AT1G01020.1
Chr1	TAIR8	CDS	8417	8464	.	-	0	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR8	exon	8417	8464	.	-	.	Parent=AT1G01020.1
Chr1	TAIR8	CDS	8236	8325	.	-	0	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR8	exon	8236	8325	.	-	.	Parent=AT1G01020.1
Chr1	TAIR8	CDS	7942	7987	.	-	0	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR8	exon	7942	7987	.	-	.	Parent=AT1G01020.1
Chr1	TAIR8	CDS	7762	7835	.	-	2	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR8	exon	7762	7835	.	-	.	Parent=AT1G01020.1
Chr1	TAIR8	CDS	7564	7649	.	-	0	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR8	exon	7564	7649	.	-	.	Parent=AT1G01020.1
Chr1	TAIR8	CDS	7384	7450	.	-	1	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR8	exon	7384	7450	.	-	.	Parent=AT1G01020.1
Chr1	TAIR8	CDS	7157	7232	.	-	0	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR8	exon	7157	7232	.	-	.	Parent=AT1G01020.1
Chr1	TAIR8	CDS	6915	7069	.	-	2	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR8	three_prime_UTR	6790	6914	.	-	.	Parent=AT1G01020.1
Chr1	TAIR8	exon	6790	7069	.	-	.	Parent=AT1G01020.1
Chr1	TAIR8	mRNA	6790	8737	.	-	.	ID=AT1G01020.2;Parent=AT1G01020;Name=AT1G01020.2;Index=1
Chr1	TAIR8	protein	7315	8666	.	-	.	ID=AT1G01020.2-Protein;Name=AT1G01020.2;Derives_from=AT1G01020.2
Chr1	TAIR8	five_prime_UTR	8667	8737	.	-	.	Parent=AT1G01020.2
Chr1	TAIR8	CDS	8571	8666	.	-	0	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR8	exon	8571	8737	.	-	.	Parent=AT1G01020.2
Chr1	TAIR8	CDS	8417	8464	.	-	0	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR8	exon	8417	8464	.	-	.	Parent=AT1G01020.2
Chr1	TAIR8	CDS	8236	8325	.	-	0	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR8	exon	8236	8325	.	-	.	Parent=AT1G01020.2
Chr1	TAIR8	CDS	7942	7987	.	-	0	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR8	exon	7942	7987	.	-	.	Parent=AT1G01020.2
Chr1	TAIR8	CDS	7762	7835	.	-	2	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR8	exon	7762	7835	.	-	.	Parent=AT1G01020.2
Chr1	TAIR8	CDS	7564	7649	.	-	0	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR8	exon	7564	7649	.	-	.	Parent=AT1G01020.2
Chr1	TAIR8	CDS	7315	7450	.	-	1	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR8	three_prime_UTR	7157	7314	.	-	.	Parent=AT1G01020.2
Chr1	TAIR8	exon	7157	7450	.	-	.	Parent=AT1G01020.2
Chr1	TAIR8	three_prime_UTR	6790	7069	.	-	.	Parent=AT1G01020.2
Chr1	TAIR8	exon	6790	7069	.	-	.	Parent=AT1G01020.2
Chr1	TAIR8	mRNA	11649	13714	.	-	.	ID=AT1G01030.1;Parent=AT1G01030;Name=AT1G01030.1;Index=1
Chr1	TAIR8	protein	11864	12940	.	-	.	ID=AT1G01030.1-Protein;Name=AT1G01030.1;Derives_from=AT1G01030.1
Chr1	TAIR8	three_prime_UTR	13335	13714	.	-	.	Parent=AT1G01030.1
Chr1	TAIR8	exon	13335	13714	.	-	.	Parent=AT1G01030.1
Chr1	TAIR8	five_prime_UTR	12941	13173	.	-	.	Parent=AT1G01030.1
Chr1	TAIR8	CDS	11864	12940	.	-	0	Parent=AT1G01030.1,AT1G01030.1-Protein;
Chr1	TAIR8	three_prime_UTR	11649	11863	.	-	.	Parent=AT1G01030.1
Chr1	TAIR8	exon	11649	13173	.	-	.	Parent=AT1G01030.1
Chr1	TAIR8	mRNA	23146	31227	.	+	.	ID=AT1G01040.1;Parent=AT1G01040;Name=AT1G01040.1;Index=1
Chr1	TAIR8	protein	23519	31079	.	+	.	ID=AT1G01040.1-Protein;Name=AT1G01040.1;Derives_from=AT1G01040.1
Chr1	TAIR8	exon	23146	24451	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	five_prime_UTR	23146	23518	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	23519	24451	.	+	0	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	24542	24655	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	24542	24655	.	+	0	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	24752	24962	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	24752	24962	.	+	0	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	25041	25435	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	25041	25435	.	+	2	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	25524	25743	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	25524	25743	.	+	0	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	25825	25997	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	25825	25997	.	+	2	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	26081	26203	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	26081	26203	.	+	0	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	26292	26452	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	26292	26452	.	+	0	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	26543	26776	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	26543	26776	.	+	1	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	26862	27012	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	26862	27012	.	+	1	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	27099	27281	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	27099	27281	.	+	0	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	27372	27533	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	27372	27533	.	+	0	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	27618	27713	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	27618	27713	.	+	0	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	27803	28431	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	27803	28431	.	+	0	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	28708	28805	.	+	.	Parent=AT1G01040.1
Chr1	TAIR8	CDS	28708	28805	.	+	1	Parent=AT1G01040.1,AT1G01040.1-Protein;
Chr1	TAIR8	exon	28890	29080	.	+	.	Parent=AT1G01040.1
