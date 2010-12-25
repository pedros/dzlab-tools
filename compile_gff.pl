#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';

use FindBin;
use lib "$FindBin::Bin/DZLab-Tools/lib";
use DZLab::Tools::GFF qw/gff_to_string gff_make_iterator/;
use DZLab::Tools::GFFStore;

use Pod::Usage;
use Getopt::Long;

my $help;
my $verbose = 0;
my $outfile;
my $memory = 0;
my $result = GetOptions (
    "verbose" => \$verbose,
    "help"    => \$help,
    "outfile|o=s" => \$outfile,
    "memory|m" => \$memory,
);
pod2usage(-verbose => 1) if (!$result || $help || !$outfile);  

unless ($outfile eq '-'){
    close STDOUT;
    open STDOUT, '>', $outfile or die "can't open $outfile for writing";
}

my $gffstore = DZLab::Tools::GFFStore->new({
        attributes => {c => 'numeric', t => 'numeric'}, 
        verbose => $verbose, 
        memory => $memory,
        #debug => 1,
    });

for (@ARGV){
    $gffstore->slurp({filename => $_});
}

my $iter = $gffstore->select_iter(<<SELECT );
    select 
    seqname,
    source,
    feature,
    start,
    end,
    (cast(sum(c) as real)/((sum(t)+sum(c)))) as score,
    strand,
    frame,
    sum(c) as c,
    sum(t) as t 
    from
    gff group by seqname, source, start, end
SELECT

while (my $row = $iter->()){
    say gff_to_string(
        [ 
        @{$row}{qw/seqname source feature start end score strand frame/},
        # make the attributes
        (join q{;}, map { $_ . "=" . $row->{$_} } qw/c t/)
        ]
    );
}

=head1 NAME

compile_gff.pl - sum the 'c' and 't' field of gff records.

=head1 SYNOPSIS

compile_gff.pl [-h] [-v] [-m] -o outfile.gff gff_file1.gff gff_file2.gff ...

=head1 DESCRIPTION

Add the 't' and 'c' attributes for gff records with identical seqname, source,
start, and end coordinates.  Score is updated to c/(t+c)

=head1 OPTIONS

    --verbose | -v          print verbose error messages
    --help    | -h          print this information
    --outfile | -o [file]   output file. use '-' to dump to screen (standard out). 
    --memory  | -m          do everything in memory (not the default). this may be faster
                            but be careful when the filesizes add up to more than RAM.

=cut

