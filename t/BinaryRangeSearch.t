#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use List::Util qw/shuffle/;
use feature 'say';
use Test::More;

use DZLab::Tools::BinaryRangeSearch;
use DZLab::Tools::GFF qw/gff_slurp_by_seq gff_to_string/;

my $gff_records = gff_slurp_by_seq {handle => \*DATA, debug => 0, sort => 1};

my $database = $gff_records->{chr1};
#my $database = [shuffle @{$gff_records->{chr1}}];


for my $gffrec (@$database){
    next unless (ref $gffrec eq 'HASH' and keys %$gffrec);
    my $target = $gffrec->{start} ;
    my $it = binary_range_search {queries => [[$target,$target]], database => $database};

    my @results;
    while (my $gff = $it->()){
        push @results, $gff;
    }
    ok(scalar @results >= 1, "at least one result found");
    #is_deeply($results[0], $gffrec);
}
#print Dumper $database;
done_testing();

__DATA__
Chr1	TAIR8	gene	3631	5899	.	+	.	ID=AT1G01010;Name=AT1G01010;Note=ANAC001 (Arabidopsis NAC domain containing protein 1),transcription factor
Chr1	TAIR8	gene	6790	8737	.	-	.	ID=AT1G01020;Name=AT1G01020;Note=ARV1
Chr1	TAIR8	gene	11649	13714	.	-	.	ID=AT1G01030;Name=AT1G01030;Note=NGA3 (NGATHA3),transcription factor
Chr1	TAIR8	gene	23146	31227	.	+	.	ID=AT1G01040;Name=AT1G01040;Note=DCL1 (DICER-LIKE1),ATP-dependent helicase,ribonuclease III
Chr1	TAIR8	gene	28500	28706	.	+	.	ID=AT1G01046;Name=AT1G01046;Note=MIR838a,miRNA
Chr1	TAIR8	gene	31170	33153	.	-	.	ID=AT1G01050;Name=AT1G01050;Note=ATPPA1 (ARABIDOPSIS THALIANA PYROPHOSPHORYLASE 1),inorganic diphosphatase,pyrophosphatase

Chr1	TAIR8	gene	33379	37840	.	-	.	ID=AT1G01060;Name=AT1G01060;Note=LHY (LATE ELONGATED HYPOCOTYL),DNA binding,transcription factor
Chr1	TAIR8	gene	38752	40944	.	-	.	ID=AT1G01070;Name=AT1G01070;Note=nodulin MtN21 family protein
Chr1	TAIR8	gene	45296	47019	.	-	.	ID=AT1G01080;Name=AT1G01080;Note=33 kDa ribonucleoprotein,chloroplast,putative,RNA-binding protein cp33,putative
Chr1	TAIR8	gene	47485	49286	.	-	.	ID=AT1G01090;Name=AT1G01090;Note=PDH-E1 ALPHA (PYRUVATE DEHYDROGENASE E1 ALPHA),pyruvate dehydrogenase (acetyl-transferring)
Chr1	TAIR8	gene	50075	51199	.	-	.	ID=AT1G01100;Name=AT1G01100;Note=60S acidic ribosomal protein P1 (RPP1A)
Chr1	TAIR8	gene	52239	54692	.	+	.	ID=AT1G01110;Name=AT1G01110;Note=IQD18 (IQ-domain 18)
Chr1	TAIR8	gene	57269	59167	.	-	.	ID=AT1G01120;Name=AT1G01120;Note=KCS1 (3-KETOACYL-COA SYNTHASE 1),acyltransferase
Chr1	TAIR8	gene	61963	63811	.	-	.	ID=AT1G01130;Name=AT1G01130
Chr1	TAIR8	gene	64166	67625	.	-	.	ID=AT1G01140;Name=AT1G01140;Note=CIPK9 (CBL-INTERACTING PROTEIN KINASE 9),kinase
Chr1	TAIR8	gene	70115	72138	.	-	.	ID=AT1G01150;Name=AT1G01150;Note=DNA binding,zinc ion binding
Chr1	TAIR8	gene	72339	74096	.	+	.	ID=AT1G01160;Name=AT1G01160;Note=GIF2 (GRF1-INTERACTING FACTOR 2)
Chr1	TAIR8	gene	73931	74737	.	-	.	ID=AT1G01170;Name=AT1G01170;Note=ozone-responsive stress-related protein,putative
Chr1	TAIR8	gene	75633	77446	.	+	.	ID=AT1G01180;Name=AT1G01180
Chr1	TAIR8	gene	78932	79032	.	-	.	ID=AT1G01183;Name=AT1G01183;Note=MIR165/MIR165A,miRNA
Chr1	TAIR8	gene	83045	84864	.	-	.	ID=AT1G01190;Name=AT1G01190;Note=CYP78A8 (cytochrome P450,family 78,subfamily A,polypeptide 8),oxygen binding
Chr1	TAIR8	gene	86515	88213	.	-	.	ID=AT1G01200;Name=AT1G01200;Note=AtRABA3 (Arabidopsis Rab GTPase homolog A3),GTP binding
Chr1	TAIR8	gene	88898	89745	.	+	.	ID=AT1G01210;Name=AT1G01210;Note=DNA-directed RNA polymerase III family protein
Chr1	TAIR8	gene	91750	95651	.	+	.	ID=AT1G01220;Name=AT1G01220;Note=GHMP kinase-related
Chr1	TAIR8	gene	95987	97407	.	+	.	ID=AT1G01225;Name=AT1G01225;Note=NC domain-containing protein-related
Chr1	TAIR8	gene	97456	99240	.	+	.	ID=AT1G01230;Name=AT1G01230;Note=ORMDL family protein
Chr1	TAIR8	gene	99894	101834	.	+	.	ID=AT1G01240;Name=AT1G01240
Chr1	TAIR8	gene	104491	105330	.	-	.	ID=AT1G01250;Name=AT1G01250;Note=AP2 domain-containing transcription factor,putative
Chr1	TAIR8	gene	109032	111609	.	+	.	ID=AT1G01260;Name=AT1G01260;Note=basic helix-loop-helix (bHLH) family protein
Chr1	TAIR8	gene	111890	111961	.	-	.	ID=AT1G01270;Name=AT1G01270;Note=pre-tRNA
Chr1	TAIR8	gene	112263	113947	.	+	.	ID=AT1G01280;Name=AT1G01280;Note=CYP703/CYP703A2 (CYTOCHROME P450,FAMILY 703,SUBFAMILY A,POLYPEPTIDE 2),oxidoreductase,acting on paired donors,with incorporation or reduction of molecular oxygen,NADH or NADPH as one donor,and incorporation of one atom of oxygen,oxygen binding
Chr1	TAIR8	gene	114286	115549	.	+	.	ID=AT1G01290;Name=AT1G01290;Note=CNX3 (COFACTOR OF NITRATE REDUCTASE AND XANTHINE DEHYDROGENASE 3),catalytic
Chr1	TAIR8	gene	116943	118764	.	+	.	ID=AT1G01300;Name=AT1G01300;Note=aspartyl protease family protein
Chr1	TAIR8	gene	119397	119997	.	+	.	ID=AT1G01305;Name=AT1G01305;Note=unknown protein
Chr1	TAIR8	gene	120154	121130	.	+	.	ID=AT1G01310;Name=AT1G01310;Note=allergen V5/Tpx-1-related family protein
Chr1	TAIR8	gene	121124	130099	.	-	.	ID=AT1G01320;Name=AT1G01320;Note=tetratricopeptide repeat (TPR)-containing protein
Chr1	TAIR8	gene	132328	135322	.	-	.	ID=AT1G01340;Name=AT1G01340;Note=ATCNGC10 (CYCLIC NUCLEOTIDE GATED CHANNEL 10),calmodulin binding,cyclic nucleotide binding,ion channel
Chr1	TAIR8	gene	136124	138162	.	+	.	ID=AT1G01350;Name=AT1G01350;Note=nucleic acid binding
Chr1	TAIR8	gene	138513	139568	.	+	.	ID=AT1G01355;Name=AT1G01355;Note=nucleic acid binding,zinc ion binding
Chr1	TAIR8	gene	141971	143183	.	+	.	ID=AT1G01360;Name=AT1G01360
Chr1	TAIR8	gene	143564	145684	.	+	.	ID=AT1G01370;Name=AT1G01370;Note=HTR12 (CENTROMERIC HISTONE H3),DNA binding
Chr1	TAIR8	gene	147153	147942	.	+	.	ID=AT1G01380;Name=AT1G01380;Note=ETC1 (ENHANCER OF TRY AND CPC 1),DNA binding,transcription factor
Chr1	TAIR8	gene	148120	149806	.	-	.	ID=AT1G01390;Name=AT1G01390;Note=UDP-glucoronosyl/UDP-glucosyl transferase family protein
Chr2	TAIR8	gene	1871	2111	.	+	.	ID=AT2G01008;Name=AT2G01008;Note=other RNA
Chr2	TAIR8	gene	3706	5513	.	+	.	ID=AT2G01010;Name=AT2G01010;Note=rRNA
Chr2	TAIR8	gene	5782	5945	.	+	.	ID=AT2G01020;Name=AT2G01020;Note=rRNA
Chr2	TAIR8	gene	6571	6672	.	+	.	ID=AT2G01021;Name=AT2G01021;Note=unknown protein
Chr2	TAIR8	gene	9648	9767	.	-	.	ID=AT2G01023;Name=AT2G01023;Note=unknown protein
Chr2	TAIR8	gene	68337	69884	.	-	.	ID=AT2G01050;Name=AT2G01050;Note=nucleic acid binding,zinc ion binding
Chr2	TAIR8	gene	73258	75389	.	-	.	ID=AT2G01060;Name=AT2G01060;Note=myb family transcription factor
Chr2	TAIR8	gene	75520	77865	.	+	.	ID=AT2G01070;Name=AT2G01070
Chr2	TAIR8	gene	77888	79334	.	+	.	ID=AT2G01080;Name=AT2G01080
Chr2	TAIR8	gene	80014	81200	.	+	.	ID=AT2G01090;Name=AT2G01090;Note=ubiquinol-cytochrome C reductase complex 7.8 kDa protein,putative,mitochondrial hinge protein,putative
Chr2	TAIR8	gene	81436	83290	.	+	.	ID=AT2G01100;Name=AT2G01100
Chr2	TAIR8	gene	83260	85229	.	-	.	ID=AT2G01110;Name=AT2G01110;Note=APG2 (ALBINO AND PALE GREEN 2)
Chr2	TAIR8	gene	85444	85608	.	+	.	ID=AT2G01111;Name=AT2G01111;Note=unknown protein
Chr2	TAIR8	gene	87317	88122	.	+	.	ID=AT2G01120;Name=AT2G01120;Note=ATORC4/ORC4 (ORIGIN RECOGNITION COMPLEX SUBUNIT 4),protein binding
Chr2	TAIR8	gene	88846	94659	.	-	.	ID=AT2G01130;Name=AT2G01130;Note=ATP binding,helicase,nucleic acid binding
Chr2	TAIR8	gene	94810	96654	.	-	.	ID=AT2G01140;Name=AT2G01140;Note=fructose-bisphosphate aldolase,putative
Chr2	TAIR8	gene	100648	101494	.	+	.	ID=AT2G01150;Name=AT2G01150;Note=RHA2B (RING-H2 FINGER PROTEIN 2B),protein binding,zinc ion binding
Chr2	TAIR8	gene	102064	102137	.	+	.	ID=AT2G01160;Name=AT2G01160;Note=pre-tRNA
Chr2	TAIR8	gene	102193	104548	.	-	.	ID=AT2G01170;Name=AT2G01170;Note=amino acid permease family protein
Chr2	TAIR8	gene	104883	106215	.	-	.	ID=AT2G01175;Name=AT2G01175;Note=unknown protein
Chr2	TAIR8	gene	106899	108802	.	-	.	ID=AT2G01180;Name=AT2G01180;Note=ATPAP1 (PHOSPHATIDIC ACID PHOSPHATASE 1),phosphatidate phosphatase
Chr2	TAIR8	gene	114974	117639	.	+	.	ID=AT2G01190;Name=AT2G01190;Note=octicosapeptide/Phox/Bem1p (PB1) domain-containing protein
Chr2	TAIR8	gene	118016	119341	.	+	.	ID=AT2G01200;Name=AT2G01200;Note=IAA32 (INDOLEACETIC ACID-INDUCED PROTEIN 32),transcription factor
Chr2	TAIR8	gene	119440	121845	.	-	.	ID=AT2G01210;Name=AT2G01210;Note=leucine-rich repeat transmembrane protein kinase,putative
Chr2	TAIR8	gene	123275	126465	.	+	.	ID=AT2G01220;Name=AT2G01220;Note=nucleotidyltransferase
Chr2	TAIR8	gene	128439	129255	.	-	.	ID=AT2G01240;Name=AT2G01240;Note=reticulon family protein (RTNLB15)
Chr2	TAIR8	gene	132696	134454	.	-	.	ID=AT2G01250;Name=AT2G01250;Note=60S ribosomal protein L7 (RPL7B)
Chr2	TAIR8	gene	135243	137805	.	-	.	ID=AT2G01260;Name=AT2G01260
Chr2	TAIR8	gene	139330	142559	.	+	.	ID=AT2G01270;Name=AT2G01270;Note=ATQSOX2 (QUIESCIN-SULFHYDRYL OXIDASE 2),thiol-disulfide exchange intermediate
Chr2	TAIR8	gene	142497	144523	.	-	.	ID=AT2G01275;Name=AT2G01275;Note=zinc finger (C3HC4-type RING finger) family protein
Chr2	TAIR8	gene	145675	149393	.	+	.	ID=AT2G01280;Name=AT2G01280;Note=MEE65 (maternal effect embryo arrest 65),RNA polymerase II transcription factor,cation:chloride symporter
