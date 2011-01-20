package SequenceLength;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw($sequence_length);

our $sequence_length = 
{ 
    arabidopsis  => {
        chr1 => 30432563,
        chr2 => 19705359,
        chr3 => 23470805,
        chr4 => 18585042,
        chr5 => 26992728,
        chrc => 154478,
        chrm => 366924,
    },
    rice => {
        chr01 => 43596771,
        chr02 => 35925388,
        chr03 => 36345490,
        chr04 => 35244269,
        chr05 => 29874162,
        chr06 => 31246789,
        chr07 => 29688601,
        chr08 => 28309179,
        chr09 => 23011239,
        chr10 => 22876596,
        chr11 => 28462103,
        chr12 => 27497214,
        chrc  => 134525,
        chrm  => 490520,
    },
    puffer => {
        1           => 22981688,
        10          => 13272281,
        11          => 11954808,
        12          => 12622881,
        13          => 13302670,
        14          => 10246949,
        15          => 7320470,
        '15_random' => 3234215,
        16          => 9031048,
        17          => 12136232,
        18          => 11077504,
        19          => 7272499,
        '1_random'  => 1180980,
        2           => 21591555,
        20          => 3798727,
        21          => 5834722,
        '21_random' => 3311323,
        '2_random'  => 2082209,
        3           => 15489435,
        4           => 9874776,
        5           => 13390619,
        6           => 7024381,
        7           => 11693588,
        8           => 10512681,
        9           => 10554956,
        mt          => 16462,
    }
};

1;
