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
use GFF::Slurp;

my $data_pos = tell DATA;

my $arr = gff_slurp(\*DATA);

is(scalar @$arr, 8, "slurp");

seek DATA, $data_pos, 0;

my $hash = gff_slurp_index(\*DATA,'sequence');

is(scalar @{$hash->{chr1}},6,"slurp index 1");
is(scalar @{$hash->{chr2}},1,"slurp index 2");

__DATA__
Chr1	TAIR8	gene	3631	5899	.	+	.	AT1G01010; ANAC001 (Arabidopsis NAC domain containing protein 1); transcription factor
Chr1	TAIR8	gene	6790	8737	.	-	.	AT1G01020; ARV1
Chr1	TAIR8	gene	11649	13714	.	-	.	AT1G01030; NGA3 (NGATHA3); transcription factor
Chr2	TAIR8	gene	23146	31227	.	+	.	AT1G01040; DCL1 (DICER-LIKE1); ATP-dependent helicase/ ribonuclease III
chr1	TAIR8	gene	28500	28706	.	+	.	AT1G01046; MIR838a; miRNA
Chr1	TAIR8	gene	31170	33153	.	-	.	AT1G01050; ATPPA1 (ARABIDOPSIS THALIANA PYROPHOSPHORYLASE 1); inorganic diphosphatase/ pyrophosphatase
Chr1	TAIR8	gene	33379	37840	.	-	.	AT1G01060; LHY (LATE ELONGATED HYPOCOTYL); DNA binding / transcription factor
.	TAIR8	gene	33379	37840	.	-	.	AT1G01060; LHY (LATE ELONGATED HYPOCOTYL); DNA binding / transcription factor
