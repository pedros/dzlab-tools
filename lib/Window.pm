#!/usr/bin/env perl
# support file for window_gff_new.pl

#==================================================================
# Window

package Window;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use autodie;
use Carp;

sub accumulate{ 
    my ($self, $gff) = @_;

    # default accumulations
    $self->counter($self->counter()+1);
    $self->score_total($self->score_total()+($gff->{score}//0));
}

sub to_gff{
    croak "abstract class: unimplemented";
}

has seqname     => (is => rw, isa => Str, required => 1);
has feature     => (is => ro, isa => Str, default  => window);
has source      => (is => ro, isa => Str, default  => dzlab);
has start       => (is => ro, isa => Int, required => 1);
has end         => (is => ro, isa => Int, required => 1);
has score_total => (is => rw, isa => Num, default  => 0);
has strand      => (is => ro, isa => Str, default  => '.');
has frame       => (is => ro, isa => Str, default  => '.');
has counter     => (is => rw, isa => Int, default  => 0);

no Moose;
__PACKAGE__->meta->make_immutable;

1;

#==================================================================
# Window::SumScore

package Window::SumScore;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use autodie;
use Carp;

extends 'Window';

# if true, reverse the count ('n' in attr) with the score... honestly not sure why this needs to be happen
# but it's in the original window_gff.pl
has 'reverse' => (is => 'ro', isa => 'Bool', default => 0);

override to_gff => sub {
    my ($self) = @_;
    my ($score, $count) = ($self->score_total,$self->counter);
    ($score,$count) = ($count,$score) if $self->reverse();
    return join "\t", 
    $self->seqname, $self->source, $self->feature, $self->start, $self->end,
    $score,
    $self->strand,
    $self->frame,
    "n=" . $count;
};

#==================================================================
# Window::FractionalMethylation

package Window::FractionalMethylation;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use autodie;
use Carp;

extends 'Window';

has t => (is => rw isa => Int default => 0);
has c => (is => rw isa => Int default => 0);

around 'accumulate' => sub {
    my ($orig,$self,$gff) = @_;
    $self->$orig($gff);
    $self->t($self->t() + $gff->{t});
    $self->c($self->c() + $gff->{c});
    # more
};

override to_gff => sub {
    my ($self) = @_;
    return join "\t", 
    $self->seqname, $self->source, $self->feature, $self->start, $self->end,
    $self->methscore(),
    $self->strand,
    $self->frame,
    sprintf("n=%d;c=%d;t=%d", $self->counter, $self->c, $self->t)
};

sub methscore{ my ($self) = @_; return $self->c() / ($self->c() + $self->t()); }

no Moose;
__PACKAGE__->meta->make_immutable;

1;

#==================================================================
# Window::Average 

package Window::AverageScore;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use autodie;
use Carp;

extends 'Window';

has score_avg => (is => rw isa => Num default => 0);
has score_std => (is => rw isa => Num default => 0);
has score_var => (is => rw isa => Num default => 0);

# is this right?
around 'accumulate' => sub {
    my ($orig,$self,$gff) = @_;
    $self->$orig($gff);

    my $previous_score_avg = $self->score_avg;

    $self->score_avg(
        $self->score_avg + ( $gff->{score} - $self->score_avg ) / $self->counter()
    );

    $self->score_std(
        $self->score_std
        + ( $gff->{score} - $previous_score_avg )
        * ( $gff->{score} - $self->score_avg )
    );

    $self->score_var(
        $self->score_std / ( $self->counter - 1 )
    ) if $self->counter > 1;
};

override to_gff => sub{
    my ($self) = @_;
    return join "\t", 
    $self->seqname, $self->source, $self->feature, $self->start, $self->end,
    $self->score_avg,
    $self->strand,
    $self->frame,
    sprintf("n=%d;var=%d;std=%d", $self->counter, $self->score_var, sqrt($self->score_std))
};

no Moose;
__PACKAGE__->meta->make_immutable;

1;

#==================================================================
# Window::Locus 

package Window::Locus;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use autodie;
use Carp;

extends 'Window';

has loci      => (is => ro isa => ArrayRef[Str] default => sub {[]});
has locus_tag => (is => ro isa => Str init_arg   => locus);

around 'accumulate' => sub {
    my ($orig,$self,$gff) = @_;
    $self->$orig($gff);
    push @{$self->loci}, $gff->{$self->locus_tag};
};

override to_gff => sub{
    my ($self) = @_;
    return join "\t", 
    $self->seqname, $self->source, $self->feature, $self->start, $self->end,
    $self->counter,
    $self->strand,
    $self->frame,
    sprintf("n=%d;%s_id=%s",
        $self->counter,
        $self->locus_tag,
        join ",", @{$self->loci}
    );
};

no Moose;
__PACKAGE__->meta->make_immutable;

1;

#==================================================================
# main


package main;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Carp;
use Smart::Comments;

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser::Attributes;

my $data_pos = tell DATA;


unless (caller()){

    ### FractionalMethylation
    {
        seek DATA, $data_pos, 0;
        my $w = Window::FractionalMethylation->new(
            seqname => 'chr1', 
            start => 1, 
            end => 50
        ); 
        my $p = GFF::Parser::Attributes->new(file => \*DATA);

        while (my $gff = $p->next()){ $w->accumulate($gff); }
        say $w->to_gff();
    }

    ### SumScore
    {
        seek DATA, $data_pos, 0;
        my $w = Window::SumScore->new( seqname => 'chr1', start => 1, end => 50, reverse => 1); 
        my $x = Window::SumScore->new( seqname => 'chr1', start => 1, end => 50); 
        my $p = GFF::Parser::Attributes->new(file => \*DATA);

        while (my $gff = $p->next()){ 
            $w->accumulate($gff); 
            $x->accumulate($gff); 
        }
        say $w->to_gff();
        say $x->to_gff();
    }

    ### AverageScore
    {
        seek DATA, $data_pos, 0;
        my $w = Window::AverageScore->new(
            seqname => 'chr1', 
            start => 1, 
            end => 50
        ); 
        my $p = GFF::Parser::Attributes->new(file => \*DATA);

        while (my $gff = $p->next()){ $w->accumulate($gff); }
        say $w->to_gff();
    }

    ### Locus
    {
        seek DATA, $data_pos, 0;
        my $w = Window::Locus->new(
            seqname => 'chr1', 
            start => 1, 
            end => 50,
            locus => 'ID'
        ); 
        my $p = GFF::Parser::Attributes->new(file => \*DATA);

        while (my $gff = $p->next()){ $w->accumulate($gff); }
        say $w->to_gff();
    }

}

__DATA__
Chr1	TAIR8	gene	3631	5899	3.14	+	.	c=1;t=3;ID=AT1G01010; ANAC001 (Arabidopsis NAC domain containing protein 1); transcription factor
Chr1	TAIR8	gene	6790	8737	12.2	-	.	c=1;t=3;ID=AT1G01020; ARV1
Chr1	TAIR8	gene	11649	13714	12.2	-	.	c=1;t=3;ID=AT1G01030; NGA3 (NGATHA3); transcription factor
Chr1	TAIR8	gene	23146	31227	3.14	+	.	c=1;t=3;ID=AT1G01040; DCL1 (DICER-LIKE1); ATP-dependent helicase/ ribonuclease III
Chr1	TAIR8	gene	28500	28706	3.14	+	.	c=1;t=3;ID=AT1G01046; MIR838a; miRNA
Chr1	TAIR8	gene	31170	33153	12.2	-	.	c=1;t=3;ID=AT1G01050; ATPPA1 (ARABIDOPSIS THALIANA PYROPHOSPHORYLASE 1); inorganic diphosphatase/ pyrophosphatase
Chr1	TAIR8	gene	33379	37840	12.2	-	.	c=1;t=3;ID=AT1G01060; LHY (LATE ELONGATED HYPOCOTYL); DNA binding / transcription factor
Chr1	TAIR8	gene	38752	40944	12.2	-	.	c=1;t=3;ID=AT1G01070; nodulin MtN21 family protein
Chr1	TAIR8	gene	45296	47019	12.2	-	.	c=1;t=3;ID="AT1G01080; 33 kDa ribonucleoprotein, chloroplast, putative / RNA-binding protein cp33, putative"
Chr1	TAIR8	gene	47485	49286	12.2	-	.	c=1;t=3;ID=AT1G01090; PDH-E1 ALPHA (PYRUVATE DEHYDROGENASE E1 ALPHA); pyruvate dehydrogenase (acetyl-transferring)
Chr1	TAIR8	gene	50075	51199	12.2	-	.	c=1;t=3;ID=AT1G01100; 60S acidic ribosomal protein P1 (RPP1A)
Chr1	TAIR8	gene	52239	54692	3.14	+	.	c=1;t=3;ID=AT1G01110; IQD18 (IQ-domain 18)
Chr1	TAIR8	gene	57269	59167	12.2	-	.	c=1;t=3;ID=AT1G01120; KCS1 (3-KETOACYL-COA SYNTHASE 1); acyltransferase
Chr1	TAIR8	gene	61963	63811	12.2	-	.	c=1;t=3;ID=AT1G01130; 0
Chr1	TAIR8	gene	64166	67625	12.2	-	.	c=1;t=3;ID=AT1G01140; CIPK9 (CBL-INTERACTING PROTEIN KINASE 9); kinase
Chr1	TAIR8	gene	70115	72138	12.2	-	.	c=1;t=3;ID=AT1G01150; DNA binding / zinc ion binding
Chr1	TAIR8	gene	72339	74096	3.14	+	.	c=1;t=3;ID=AT1G01160; GIF2 (GRF1-INTERACTING FACTOR 2)
Chr1	TAIR8	gene	73931	74737	12.2	-	.	c=1;t=3;ID="AT1G01170; ozone-responsive stress-related protein, putative"
Chr1	TAIR8	gene	75633	77446	3.14	+	.	c=1;t=3;ID=AT1G01180; 0
Chr1	TAIR8	gene	78932	79032	12.2	-	.	c=1;t=3;ID=AT1G01183; MIR165/MIR165A; miRNA
Chr1	TAIR8	gene	83045	84864	12.2	-	.	c=1;t=3;ID="AT1G01190; CYP78A8 (cytochrome P450, family 78, subfamily A, polypeptide 8); oxygen binding"
Chr1	TAIR8	gene	86515	88213	12.2	-	.	c=1;t=3;ID=AT1G01200; AtRABA3 (Arabidopsis Rab GTPase homolog A3); GTP binding
Chr1	TAIR8	gene	88898	89745	3.14	+	.	c=1;t=3;ID=AT1G01210; DNA-directed RNA polymerase III family protein
Chr1	TAIR8	gene	91750	95651	3.14	+	.	c=1;t=3;ID=AT1G01220; GHMP kinase-related
Chr1	TAIR8	gene	95987	97407	3.14	+	.	c=1;t=3;ID=AT1G01225; NC domain-containing protein-related
Chr1	TAIR8	gene	97456	99240	3.14	+	.	c=1;t=3;ID=AT1G01230; ORMDL family protein
Chr1	TAIR8	gene	99894	101834	3.14	+	.	c=1;t=3;ID=AT1G01240; 0
Chr1	TAIR8	gene	104491	105330	12.2	-	.	c=1;t=3;ID="AT1G01250; AP2 domain-containing transcription factor, putative"
Chr1	TAIR8	gene	109032	111609	3.14	+	.	c=1;t=3;ID=AT1G01260; basic helix-loop-helix (bHLH) family protein
Chr1	TAIR8	gene	111890	111961	12.2	-	.	c=1;t=3;ID=AT1G01270; pre-tRNA
Chr1	TAIR8	gene	112263	113947	3.14	+	.	c=1;t=3;ID="AT1G01280; CYP703/CYP703A2 (CYTOCHROME P450, FAMILY 703, SUBFAMILY A, POLYPEPTIDE 2); oxidoreductase, acting on paired donors, with incorporation or reduction of molecular oxygen, NADH or NADPH as one donor, and incorporation of one atom of oxygen / oxygen binding"
Chr1	TAIR8	gene	114286	115549	3.14	+	.	c=1;t=3;ID=AT1G01290; CNX3 (COFACTOR OF NITRATE REDUCTASE AND XANTHINE DEHYDROGENASE 3); catalytic
Chr1	TAIR8	gene	116943	118764	3.14	+	.	c=1;t=3;ID=AT1G01300; aspartyl protease family protein
Chr1	TAIR8	gene	119397	119997	3.14	+	.	c=1;t=3;ID=AT1G01305; unknown protein
Chr1	TAIR8	gene	120154	121130	3.14	+	.	c=1;t=3;ID=AT1G01310; allergen V5/Tpx-1-related family protein
Chr1	TAIR8	gene	121124	130099	12.2	-	.	c=1;t=3;ID=AT1G01320; tetratricopeptide repeat (TPR)-containing protein
Chr1	TAIR8	gene	132328	135322	12.2	-	.	c=1;t=3;ID=AT1G01340; ATCNGC10 (CYCLIC NUCLEOTIDE GATED CHANNEL 10); calmodulin binding / cyclic nucleotide binding / ion channel
Chr1	TAIR8	gene	136124	138162	3.14	+	.	c=1;t=3;ID=AT1G01350; nucleic acid binding
Chr1	TAIR8	gene	138513	139568	3.14	+	.	c=1;t=3;ID=AT1G01355; nucleic acid binding / zinc ion binding
Chr1	TAIR8	gene	141971	143183	3.14	+	.	c=1;t=3;ID=AT1G01360; 0
Chr1	TAIR8	gene	143564	145684	3.14	+	.	c=1;t=3;ID=AT1G01370; HTR12 (CENTROMERIC HISTONE H3); DNA binding
Chr1	TAIR8	gene	147153	147942	3.14	+	.	c=1;t=3;ID=AT1G01380; ETC1 (ENHANCER OF TRY AND CPC 1); DNA binding / transcription factor
Chr1	TAIR8	gene	148120	149806	12.2	-	.	c=1;t=3;ID=AT1G01390; UDP-glucoronosyl/UDP-glucosyl transferase family protein
Chr1	TAIR8	gene	150689	152210	12.2	-	.	c=1;t=3;ID=AT1G01400; unknown protein
Chr1	TAIR8	gene	153113	154198	3.14	+	.	c=1;t=3;ID=AT1G01410; APUM22 (ARABIDOPSIS PUMILIO 22); RNA binding / binding
Chr1	TAIR8	gene	154492	156011	12.2	-	.	c=1;t=3;ID=AT1G01420; UDP-glucoronosyl/UDP-glucosyl transferase family protein
Chr1	TAIR8	gene	156801	158655	12.2	-	.	c=1;t=3;ID=AT1G01430; 0
Chr1	TAIR8	gene	159856	162572	12.2	-	.	c=1;t=3;ID=AT1G01440; extra-large G-protein-related
Chr1	TAIR8	gene	163419	166239	3.14	+	.	c=1;t=3;ID=AT1G01448; other RNA
Chr1	TAIR8	gene	164105	165517	12.2	-	.	c=1;t=3;ID=AT1G01450; protein kinase-related
Chr1	TAIR8	gene	166589	167842	12.2	-	.	c=1;t=3;ID=AT1G01453; 0
Chr1	TAIR8	gene	168723	171165	3.14	+	.	c=1;t=3;ID=AT1G01460; ATPIPK11; 1-phosphatidylinositol-4-phosphate 5-kinase
Chr1	TAIR8	gene	172146	172948	12.2	-	.	c=1;t=3;ID=AT1G01470; LEA14 (LATE EMBRYOGENESIS ABUNDANT 14)
Chr1	TAIR8	gene	173251	173466	3.14	+	.	c=1;t=3;ID=AT1G01471; unknown protein
Chr1	TAIR8	gene	175782	178400	3.14	+	.	c=1;t=3;ID=AT1G01480; ACS2 (1-Amino-cyclopropane-1-carboxylate synthase 2)
Chr1	TAIR8	gene	180059	182358	12.2	-	.	c=1;t=3;ID=AT1G01490; heavy-metal-associated domain-containing protein
Chr1	TAIR8	gene	185133	186923	3.14	+	.	c=1;t=3;ID=AT1G01500; 0
Chr1	TAIR8	gene	187211	190056	3.14	+	.	c=1;t=3;ID=AT1G01510; AN (ANGUSTIFOLIA)
Chr1	TAIR8	gene	190596	192139	3.14	+	.	c=1;t=3;ID=AT1G01520; myb family transcription factor
Chr1	TAIR8	gene	192640	193670	12.2	-	.	c=1;t=3;ID=AT1G01530; AGL28 (AGAMOUS-LIKE 28); DNA binding / transcription factor
Chr1	TAIR8	gene	195780	198684	3.14	+	.	c=1;t=3;ID=AT1G01540; protein kinase family protein
Chr1	TAIR8	gene	199639	201775	3.14	+	.	c=1;t=3;ID=AT1G01550; BPS1 (BYPASS 1)
Chr1	TAIR8	gene	202136	204335	3.14	+	.	c=1;t=3;ID=AT1G01560; ATMPK11 (Arabidopsis thaliana MAP kinase 11); MAP kinase/ kinase
Chr1	TAIR8	gene	205176	207435	3.14	+	.	c=1;t=3;ID=AT1G01570; fringe-related protein
Chr1	TAIR8	gene	209395	213041	3.14	+	.	c=1;t=3;ID=AT1G01580; FRO2 (FERRIC REDUCTION OXIDASE 2); ferric-chelate reductase
Chr1	TAIR8	gene	214229	217304	3.14	+	.	c=1;t=3;ID=AT1G01590; FRO1 (FERRIC REDUCTION OXIDASE 1); ferric-chelate reductase
Chr1	TAIR8	gene	218994	221286	3.14	+	.	c=1;t=3;ID="AT1G01600; CYP86A4 (cytochrome P450, family 86, subfamily A, polypeptide 4); oxygen binding"
Chr1	TAIR8	gene	221691	224340	12.2	-	.	c=1;t=3;ID=AT1G01610; ATGPAT4/GPAT4 (GLYCEROL-3-PHOSPHATE ACYLTRANSFERASE 4); 1-acylglycerol-3-phosphate O-acyltransferase/ acyltransferase
Chr1	TAIR8	gene	225665	227302	12.2	-	.	c=1;t=3;ID=AT1G01620; PIP1C (PLASMA MEMBRANE INTRINSIC PROTEIN 1;3); water channel
Chr1	TAIR8	gene	229013	230917	3.14	+	.	c=1;t=3;ID="AT1G01630; SEC14 cytosolic factor, putative / phosphoglyceride transfer protein, putative"
Chr1	TAIR8	gene	230994	232508	12.2	-	.	c=1;t=3;ID=AT1G01640; speckle-type POZ protein-related
Chr1	TAIR8	gene	232841	237817	12.2	-	.	c=1;t=3;ID=AT1G01650; peptidase
Chr1	TAIR8	gene	240057	242608	12.2	-	.	c=1;t=3;ID=AT1G01660; U-box domain-containing protein
Chr1	TAIR8	gene	242843	245988	12.2	-	.	c=1;t=3;ID=AT1G01670; U-box domain-containing protein
Chr1	TAIR8	gene	246411	248367	12.2	-	.	c=1;t=3;ID=AT1G01680; U-box domain-containing protein
Chr1	TAIR8	gene	249141	252462	3.14	+	.	c=1;t=3;ID=AT1G01690; 0
Chr1	TAIR8	gene	252947	254495	3.14	+	.	c=1;t=3;ID=AT1G01695; 0
Chr1	TAIR8	gene	259495	261474	12.2	-	.	c=1;t=3;ID=AT1G01700; ATROPGEF2/ROPGEF2 (KINASE PARTNER PROTEIN-LIKE); Rho guanyl-nucleotide exchange factor/ / protein binding
Chr1	TAIR8	gene	262828	267771	3.14	+	.	c=1;t=3;ID=AT1G01710; acyl-CoA thioesterase family protein
Chr1	TAIR8	gene	268330	269819	3.14	+	.	c=1;t=3;ID=AT1G01720; ATAF1 (Arabidopsis NAC domain containing protein 2); transcription factor
Chr1	TAIR8	gene	269792	270775	12.2	-	.	c=1;t=3;ID=AT1G01725; 0
Chr1	TAIR8	gene	270956	272061	3.14	+	.	c=1;t=3;ID=AT1G01730; 0
Chr1	TAIR8	gene	272111	274239	12.2	-	.	c=1;t=3;ID=AT1G01740; protein kinase family protein
Chr1	TAIR8	gene	275366	276310	3.14	+	.	c=1;t=3;ID="AT1G01750; actin-depolymerizing factor, putative"
Chr1	TAIR8	gene	276266	278448	12.2	-	.	c=1;t=3;ID=AT1G01760; RNA binding / adenosine deaminase
Chr1	TAIR8	gene	278615	282891	3.14	+	.	c=1;t=3;ID=AT1G01770; 0
Chr1	TAIR8	gene	282761	284245	3.14	+	.	c=1;t=3;ID=AT1G01780; LIM domain-containing protein
Chr1	TAIR8	gene	284781	291094	3.14	+	.	c=1;t=3;ID=AT1G01790; KEA1 (K EFFLUX ANTIPORTER 1); potassium:hydrogen antiporter
Chr1	TAIR8	gene	293342	295040	3.14	+	.	c=1;t=3;ID=AT1G01800; short-chain dehydrogenase/reductase (SDR) family protein
Chr1	TAIR8	gene	295221	295859	3.14	+	.	c=1;t=3;ID=AT1G01810; unknown protein
Chr1	TAIR8	gene	296001	298120	12.2	-	.	c=1;t=3;ID=AT1G01820; PEX11C
Chr1	TAIR8	gene	298535	302315	12.2	-	.	c=1;t=3;ID=AT1G01830; armadillo/beta-catenin repeat family protein
Chr1	TAIR8	gene	303537	304358	3.14	+	.	c=1;t=3;ID=AT1G01840; 0
Chr1	TAIR8	gene	304133	306287	12.2	-	.	c=1;t=3;ID=AT1G01860; PFC1 (PALEFACE 1)
Chr1	TAIR8	gene	306384	306456	3.14	+	.	c=1;t=3;ID=AT1G01870; pre-tRNA
Chr1	TAIR8	gene	306558	308991	12.2	-	.	c=1;t=3;ID="AT1G01880; DNA repair protein, putative"
Chr1	TAIR8	gene	309275	309347	12.2	-	.	c=1;t=3;ID=AT1G01890; pre-tRNA
Chr1	TAIR8	gene	310316	313130	3.14	+	.	c=1;t=3;ID=AT1G01900; ATSBT1.1; subtilase
Chr1	TAIR8	gene	313101	315902	12.2	-	.	c=1;t=3;ID="AT1G01910; anion-transporting ATPase, putative"
Chr1	TAIR8	gene	316128	319650	3.14	+	.	c=1;t=3;ID=AT1G01920; SET domain-containing protein
