#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use List::Util qw/sum/;
use Log::Log4perl qw/get_logger/;
use File::Spec::Functions;
use File::Basename;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use FindBin;
use lib "$FindBin::Bin/lib";
use Fasta qw/bisulfite_convert/;
use Launch;


pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
unless $opt_output_directory && $opt_reference_a && $opt_reference_b && $opt_raw && $opt_ecotype_a && $opt_ecotype_b && scalar %opt_splice;


my $conf=qq/
    log4perl.logger          = DEBUG, Print
    log4perl.logger.Script   = DEBUG
    log4perl.logger.PipeLine = DEBUG
    log4perl.logger.Module   = DEBUG

    log4perl.appender.Print                          = Log::Log4perl::Appender::Screen
    log4perl.appender.Print.layout                   = PatternLayout
    log4perl.appender.Print.layout.ConversionPattern = %d{HH:mm:ss} %p> (%L) %M - %m%n

    log4perl.appender.File                          = Log::Log4perl::Appender::File
    log4perl.appender.File.filename                 = $opt_raw.log
    log4perl.appender.File.layout                   = PatternLayout
    log4perl.appender.File.layout.ConversionPattern = %d{HH:mm:ss} %p> (%L) %M - %m%n
/;
Log::Log4perl::init( \$conf );

my $logger = get_logger("PipeLine");

#################################################################################

sub _gen_files {
    my ($base, $ext, @groups) = @_;

    my %split_files;
    for my $group (@groups) {
        $split_files{$group} = "$base-$group.$ext";
    }
    return \%split_files;
}

#################################################################################

# already rc'd and bs'd
my ($trim5, $trim3) = ($opt_splice{start} - 1, $opt_read_length - $opt_splice{end});
#my @chromosomes;

#$logger->info("Chromosomes: " . join ",", @chromosomes);

#######################################################################
# Handle options

$logger->info("raw file: $opt_raw");
$logger->info("reference A: $opt_reference_a");
$logger->info("reference B: $opt_reference_b");
$logger->info("ecotype A: $opt_ecotype_a");
$logger->info("ecotype B: $opt_ecotype_b");
$logger->info("splice: " . join ',', @opt_splice{qw/start end/});

mkdir $opt_output_directory;
if (! -d $opt_output_directory){
    $logger->logdie("can't create $opt_output_directory");
}

my $basename = $opt_basename;
if (! defined $basename){
    $basename = $opt_raw;
    $basename =~ s/\.\w+$//;
    $basename = basename $basename;
}
$basename = catfile($opt_output_directory,$basename);

$logger->info("basename: $basename");

#######################################################################
# BSRC genomes, build bowtie indices

my $bsrc_reference_a = $opt_reference_a . ".bsrc";
my $bsrc_reference_b = $opt_reference_b . ".bsrc";

launch("perl -S fasta_bsrc.pl $opt_reference_a > $bsrc_reference_a", 
    expected => $bsrc_reference_a, force => $opt_force >= 2 );
launch("perl -S fasta_bsrc.pl $opt_reference_b > $bsrc_reference_b", 
    expected => $bsrc_reference_b, force => $opt_force >= 2 );

for my $ref ($bsrc_reference_a,$bsrc_reference_b) {
    launch("bowtie-build $ref $ref", 
        expected => ("$ref.1.ebwt"), force => $opt_force >= 2);
}

# raw: fastq -> fasta,  c2t
$logger->info("converting fastq->fasta, then converting");
my $rawfas = "$opt_raw.fasta";
my $rawc2t = "$opt_raw.fasta.c2t";
launch("perl -S fq_all2std.pl fq2fa $opt_raw > $rawfas", expected => $rawfas, force => $opt_force>=2);
launch("perl -S convert.pl c2t $rawfas > $rawc2t", expected => $rawc2t, force => $opt_force>=2);


#######################################################################
# Run Bowtie

# run bowtie of rawc2t vs output a and b
$logger->info("running $opt_raw against eco a and b bowtie");
my $basename_a = "$basename-vs-$opt_ecotype_a";
my $basename_b = "$basename-vs-$opt_ecotype_b";
my $bowtie_a = "$basename_a.bowtie";
my $bowtie_b = "$basename_b.bowtie";
launch("bowtie $bsrc_reference_a -f -B 1 -v $opt_bowtie_mismatches --norc --best -5 $trim5 -3 $trim3 $rawc2t $bowtie_a",
    expected => $bowtie_a, force => $opt_force);
launch("bowtie $bsrc_reference_b -f -B 1 -v $opt_bowtie_mismatches --norc --best -5 $trim5 -3 $trim3 $rawc2t $bowtie_b",
    expected => $bowtie_b, force => $opt_force);

#######################################################################
# parse bowtie

$logger->info("parse bowtie -> eland");
my $eland_a = "$basename_a.eland";
my $eland_b = "$basename_b.eland";
launch("perl -S parse_bowtie.pl -u $rawfas -s @opt_splice{qw/start end/} $bowtie_a -o $eland_a",
    expected => $eland_a, force => $opt_force);
launch("perl -S parse_bowtie.pl -u $rawfas -s @opt_splice{qw/start end/} $bowtie_b -o $eland_b",
    expected => $eland_b, force => $opt_force);

#######################################################################
# Split on mismatches

$logger->info("sort the eland files");
my $eland_sorted_a = "$eland_a.sorted";
my $eland_sorted_b = "$eland_b.sorted";
launch("sort -k 1,1 -k 4,4 -S 15% $eland_a -o $eland_sorted_a",expected => $eland_sorted_a, force => $opt_force);
launch("sort -k 1,1 -k 4,4 -S 15% $eland_b -o $eland_sorted_b",expected => $eland_sorted_b, force => $opt_force);

$logger->info("split_on_mismatch_2.pl");
my $eland_filtered_a = "$eland_a.filtered";
my $eland_filtered_b = "$eland_b.filtered";
launch("perl -S split_on_mismatches_2.pl -c -ia $eland_a -ib $eland_b -oa $eland_filtered_a -ob $eland_filtered_b",
    expected => [ $eland_filtered_a, $eland_filtered_b]);


#######################################################################
# Count and calc ratios and stuff

my %counts = ($opt_ecotype_a => {}, $opt_ecotype_b => {});
my @chromosomes = qw/ chr3 chr1 chr4 chrc chr2 chrm chr5/;
for my $c (@chromosomes) {
    for my $mm (0..$opt_bowtie_mismatches, 'total') {
        $counts{$opt_ecotype_a}{$c}{$mm} = 0;
        $counts{$opt_ecotype_b}{$c}{$mm} = 0;
    }
}
$logger->info(Dumper \%counts);

for my $eco ([$opt_ecotype_a, $eland_filtered_a], [$opt_ecotype_b, $eland_filtered_b]){
    my ($ecotype, $file) = @$eco;
    open my $fh, '<', $file or die "can't open $file";
    while (defined(my $line = <$fh>)){
        my @split = split /\t/, $line;
        $split[3] =~ /(?:RC_)?(chr.*?):.*?(\d)$/i;
        my $c = lc $1;
        $counts{$ecotype}{$c}{$2}++;
        $counts{$ecotype}{$c}{total}++;
    }
    close $fh;
}

$logger->info(Dumper \%counts);

my $ratio_log = $basename . ".log";
open my $ratiofh, '>', $ratio_log;

say $ratiofh "$opt_ecotype_a / $opt_ecotype_b";

for my $mm (0 .. $opt_bowtie_mismatches, 'total') {
    say $ratiofh "$mm mismatches";
    for my $c (sort @chromosomes) {
        say $ratiofh "$c: \t $counts{$opt_ecotype_a}{$c}{$mm} / $counts{$opt_ecotype_b}{$c}{$mm} = "  . 
        ( $counts{$opt_ecotype_b}{$c}{$mm} > 0 ?
            ($counts{$opt_ecotype_a}{$c}{$mm} / $counts{$opt_ecotype_b}{$c}{$mm})
            : "inf" );
    }
    my $total_a = sum map { $counts{$opt_ecotype_a}{$_}{$mm} } @chromosomes;
    my $total_b = sum map { $counts{$opt_ecotype_b}{$_}{$mm} } @chromosomes;
    say $ratiofh "total: \t $total_a / $total_b = " 
    . ( $total_b > 0 ?  $total_a / $total_b  : 'inf');
}


#$logger->info("correlating");
#my $eland_correlated_a = "$eland_sorted_a.correlated";
#my $eland_correlated_b = "$eland_sorted_b.correlated";
#run_cmd("perl -S correlatePairedEnds.pl -l $eland_sorted_a -f $opt_original_reference_a -s $opt_read_length -1 1 -2 0 -a 1 -o $eland_correlated_a") 
#    unless (file_exists("$eland_correlated_a"));
#run_cmd("perl -S correlatePairedEnds.pl -l $eland_sorted_b --reference $opt_original_reference_b -s $opt_read_length -1 1 -2 0 -a 1 -o $eland_correlated_b")
#    unless (file_exists("$eland_correlated_b"));

# -d, -t, does not matter when no -r. no -m (max hits)

#$logger->info("split correlated");
#my @correlated_split_a => _gen_files ($eland_correlated_a, 'gff',  @chromosomes);
#my @correlated_split_b => _gen_files ($eland_correlated_b, 'gff',  @chromosomes);
#$logger->info("correlated_split_a == " . join ",", @correlated_split_a);
#$logger->info("correlated_split_b == " . join ",", @correlated_split_b);
#
# <output dir>/single-c/basename-sequence.single-c.gff
# freq  => _gen_files (File::Spec->catfile ($outdir, 'single-c', $base_name), 'single-c.gff', @groups),

# <output dir>/single-c/basename-sequence.single-c-CG.gff, 
# <output dir>/single-c/basename-sequence.single-c-CHG.gff, etc.
#cont  => [map { _gen_files (File::Spec->catfile ($outdir, 'single-c', $base_name),
#        "single-c-$_.gff", @groups) } @contexts],

# <output dir>/windows/basename-sequence.w50-CG.gff, 
# <output dir>/windows/basename-sequence.w50-CHG.gff, etc.
#wcont => [map { _gen_files (File::Spec->catfile ($out_dir, 'windows',
#            $base_name), "w${window_size}-$_.gff", @groups) } @contexts],


=head1 NAME

ratio.pl - Your program here

=head1 SYNOPSIS

Usage examples:

 ratio.pl -r raw.fastq -ea Col -b Ler -ra genome-a.fasta -rb genome-b.fasta -l 100 -m 2 -s 1 50 -o outdir -b basename

=head1 REQUIRED ARGUMENTS

=over

=back

=head1 OPTIONS

=over

=item  -r <file> | --raw <file>

FastQ format reads

=for Euclid
    file.type:        readable

=item  -ra <fasta> | --reference-a <fasta>

Genome reference for A.

=for Euclid 
    fasta.type: readable

=item  -rb <fasta> | --reference-b <fasta>

Genome reference for B.

=for Euclid 
    fasta.type: readable

=item  -ea <eco> | --ecotype-a <eco>

Ecotype A label.

=item  -eb <eco> | --ecotype-b <eco>

Ecotype B label.

=item  -s <start> <end> | --splice <start> <end>

=for Euclid
    start.type:     int
    end.type:     int

=item  -l <len> | --read-length <len>

Number of bp in each read.

=for Euclid
    len.default: 100

=item  -m <num> | --bowtie-mismatches <num>

Number of mismatches to allow in bowtie.

=for Euclid
    num.default:     2
    num.type:        int, num >= 0 && num <= 3

=item -o <dir> | --output-directory <dir>

Output Directory.

=item  -b <name> | --basename <name>

Prefix for the file names.

=item --help | -h

=item  --force <level>

Level of forcefulness in doing jobs.  1 = Redo all run-specifics.  2 = Redo bowtie-build as well.

=for Euclid
    level.default:     0
    level.type:        int, level >= 0

=back

=cut

