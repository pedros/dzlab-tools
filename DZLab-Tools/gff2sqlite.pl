#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';

use FindBin;
use lib "$FindBin::Bin/lib";
use DZLab::Tools::GFF qw/gff_make_iterator/;
use DZLab::Tools::GFFStore;

use Pod::Usage;
use Getopt::Long::Descriptive;

my $output;
my $help;

my ($opt, $usage) = describe_options(
    "Usage: %c %o ...",
    [],
    [ "output|o=s"     , "output sqlite database file"                 , { required => 1 }]    , 
    [ "commitsize|s=i" , "number of records to commit per transaction" , { default => 10000 }] , 
    [ 'bench|b'        , "benchmark: only read into memory" ]          , 
    [ 'debug|g'        , "Be debugful" ]                               , 
    [ 'quiet|q'        , "Be quiet" ]                                  , 
    [ 'help|h'         , "print this message and exit" ]               , 
);
                
print($usage->text), exit if ($opt->help);

if (!@ARGV) { pod2usage(-verbose => 1) }

my $gffstore = DZLab::Tools::GFFStore->new({
        dbname     => $opt->output,
        verbose    => ! $opt->quiet,
        commitsize => $opt->commitsize,
    });

$gffstore->slurp({handle => \*ARGV});

=head1 NAME

gff2sqlite.pl - convert gff file to sqlite.pl

=head1 SYNOPSIS

gff2sqlite.pl -o output.sqlite file1.gff [file2.gff ...]

=head1 DESCRIPTION

Slurp up given files into output sqlite database via DZLab::Tools::GFFStore module

=head1 OPTIONS

    --help         print this message
    --output, -o   output file

=cut

