#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';

use FindBin;
use lib "$FindBin::Bin/../DZLab-Tools/lib";
use DZLab::Tools::GFF qw/gff_make_iterator/;
use DZLab::Tools::GFFStore;

use Pod::Usage;
use Getopt::Long;

my $output;
my $help;

my $result = GetOptions (
    "output|o=s" => \$output,
    "help"    => \$help,
);
pod2usage(-verbose => 1) if (!$result || $help || !$output);  

if (!@ARGV) { pod2usage(-verbose => 1) }

my $gffstore = DZLab::Tools::GFFStore->new({
        dbname => $output,
        verbose => 1,
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

