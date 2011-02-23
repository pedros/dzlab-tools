package GFF::SplitSort;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(gff_splitsort);

use GFF::Sort;
use GFF::Split;

sub gff_splitsort{
    my (%opt) = @_;
    my ($file, $sequence,$feature,$colref,$tmp) = @opt{qw/file sequence feature cols tmpdir/};
    my %split = gff_split(file => $file,sequence => $sequence,feature => $feature, tmpdir => $tmp);
    for my $split_file (values %split){
        gff_sort(file => $split_file, overwrite => 1, cols => $colref, manual => 0);
    }
    return %split;
}

1;

