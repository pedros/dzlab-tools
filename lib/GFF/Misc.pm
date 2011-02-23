package GFF::Misc;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(colname2num);

my %cols = (
    sequence         => 1,
    source           => 2,
    feature          => 3,
    start            => 4,
    end              => 5,
    score            => 6,
    strand           => 7,
    frame            => 8,
    attribute_string => 9,
);

sub colname2num{
    if (exists $cols{$_[0]}){
        return $cols{$_[0]};
    } else{
        confess "no such column";
    }
}

1;

