#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

my $eland1 = $ARGV[0];
my $eland2 = $ARGV[1];
my $correl = $ARGV[2];

my $tissue = $ARGV[3];
my $batch  = 'batch-' . $ARGV[4];

my (%ecounts1, %ecounts2) = ();

open my $EL1, '<', $eland1 or croak "Can't open $eland1";
open my $EL2, '<', $eland2 or croak "Can't open $eland2";
ELAND_LINE:
while (1) {

    my $eline1 = <$EL1>;
    my $eline2 = <$EL2>;

    last ELAND_LINE if !defined $eline1 or !defined $eline2;

    my @efields1 = split /\t/, $eline1;
    my @efields2 = split /\t/, $eline2;

    ++$ecounts1{ALL};
    ++$ecounts2{ALL};

    ++$ecounts1{NM} if $efields1[2] =~ m/NM/;
    ++$ecounts2{NM} if $efields2[2] =~ m/NM/;

    ++$ecounts1{U0} if $efields1[2] =~ m/^1:[09]:[09]/;
    ++$ecounts1{U1} if $efields1[2] =~ m/^[09]:1:[09]/;
    ++$ecounts1{U2} if $efields1[2] =~ m/^[09]:[09]:1/;

    ++$ecounts2{U0} if $efields2[2] =~ m/^1:[09]:[09]/;
    ++$ecounts2{U1} if $efields2[2] =~ m/^[09]:1:[09]/;
    ++$ecounts2{U2} if $efields2[2] =~ m/^[09]:[09]:1/;

    ++$ecounts1{RT} if $efields1[2] !~ m/NM/ and $efields1[3] =~ m/,/;
    ++$ecounts2{RT} if $efields2[2] !~ m/NM/ and $efields2[3] =~ m/,/;
}
close $EL1 or croak "Can't close $eland1";
close $EL2 or croak "Can't close $eland2";

$ecounts1{UT} = $ecounts1{U0} + $ecounts1{U1} + $ecounts1{U2};
$ecounts2{UT} = $ecounts2{U0} + $ecounts2{U1} + $ecounts2{U2};

my %ccounts = ();

if ($correl) {
    open my $CORR, '<', $correl or croak "Can't open $correl";
  CORREL_LINE:
    while (<$CORR>) {
        my @cfields = split /\t/, $_;

        ++$ccounts{ALL};

        ++$ccounts{NM} if $cfields[5]  < 1;
        ++$ccounts{UT} if $cfields[5] == 1;
        ++$ccounts{RT} if $cfields[5]  > 1;

        $ccounts{RR} += 1/2 if $cfields[1] =~ m{U\/R|R\/U} and $cfields[5] == 1;
        ++$ccounts{RR} if $cfields[1] =~ m{R\/R} and $cfields[5] == 1;
        ++$ccounts{UU} if $cfields[1] =~ m{U\/U} and $cfields[5] == 1;
    }
}

open my $LOG, '>>', "${correl}-all-alignment.table" or croak "Can't write to all.alignmentstats.";

print $LOG join ("\t",
                 'TISSUE/BATCH',
                 '/1_TOT',
                 '/1_U0',
                 '/1_UT',
                 '/1_U0/UT',
                 '/2_TOT',
                 '/2_U0',
                 '/2_UT',
                 '/2_U0/UT',
                 'REPEAT_RESOLV',
                 'TOTAL_MATCHES',
                 'RR/TM'
             ), "\n" if 1;

print $LOG join ("\t",
                 join (q{|}, $tissue, $batch),
                 $ecounts1{ALL},
                 $ecounts1{U0},
                 $ecounts1{UT},
                 ($ecounts1{U0}/$ecounts1{UT}),
                 $ecounts2{ALL},
                 $ecounts2{U0},
                 $ecounts2{UT},
                 ($ecounts2{U0}/$ecounts2{UT}),
                 $ccounts{RR},
                 $ccounts{UT},
                 ($ccounts{RR} / $ccounts{UT}),
             ), "\n";


print "================================================================================\n";
print "$eland1\n";
print "================================================================================\n";
print "No. total reads:\t$ecounts1{ALL}\n";
print "No. non-matched reads:\t$ecounts1{NM}\n";
print "No. 0-mm matched reads:\t$ecounts1{U0}\n";
print "No. 1-mm matched reads:\t$ecounts1{U1}\n";
print "No. 2-mm matched reads:\t$ecounts1{U2}\n";
print "No. matched reads:\t$ecounts1{UT}\n";
print "No. repeat reads:\t$ecounts1{RT}\n";
print "Ratio of matched reads in total reads:\t", $ecounts1{UT} / $ecounts1{ALL}, "\n";
print "Ratio of 0-mm reads in total reads:\t", $ecounts1{U0} / $ecounts1{ALL}, "\n";
print "Ratio of 0-mm reads in matched reads:\t", $ecounts1{U0} / $ecounts1{UT}, "\n";
print "Ratio of non-matched reads in total reads:\t", $ecounts1{NM} / $ecounts1{ALL}, "\n";
print "Ratio of repeats reads in total reads:\t", $ecounts1{RT} / $ecounts1{ALL}, "\n";
print "================================================================================\n";

print "================================================================================\n";
print "$eland2\n";
print "================================================================================\n";
print "No. total reads:\t$ecounts2{ALL}\n";
print "No. non-matched reads:\t$ecounts2{NM}\n";
print "No. 0-mm matched reads:\t$ecounts2{U0}\n";
print "No. 2-mm matched reads:\t$ecounts2{U1}\n";
print "No. 2-mm matched reads:\t$ecounts2{U2}\n";
print "No. matched reads:\t$ecounts2{UT}\n";
print "No. repeat reads:\t$ecounts2{RT}\n";
print "Ratio of matched reads in total reads:\t", $ecounts2{UT} / $ecounts2{ALL}, "\n";
print "Ratio of 0-mm reads in total reads:\t", $ecounts2{U0} / $ecounts2{ALL}, "\n";
print "Ratio of 0-mm reads in matched reads:\t", $ecounts2{U0} / $ecounts2{UT}, "\n";
print "Ratio of non-matched reads in total reads:\t", $ecounts2{NM} / $ecounts2{ALL}, "\n";
print "Ratio of repeats reads in total reads:\t", $ecounts2{RT} / $ecounts2{ALL}, "\n";
print "================================================================================\n";

if ($correl) {
    print "================================================================================\n";
    print "$correl\n";
    print "================================================================================\n";
    print "No. total reads:\t$ccounts{ALL}\n";
    print "No. non-matched reads:\t$ccounts{NM}\n";
    print "No. matched reads:\t$ccounts{UT}\n";
    print "No. preserved unique reads:\t$ccounts{UU}\n";
    print "No. repeat reads:\t$ccounts{RT}\n"; #
    print "No. resolved repeats:\t$ccounts{RR}\n";
    print "Ratio of preserved unique reads in total matched reads:\t", $ccounts{UU} / $ccounts{UT}, "\n";
    print "Ratio of matched reads in total reads:\t", $ccounts{UT} / $ccounts{ALL}, "\n";
    print "Ratio of non-matched reads in total reads:\t", $ccounts{NM} / $ccounts{ALL}, "\n";
    print "Ratio of repeat reads in total reads:\t", $ccounts{RT} / $ccounts{ALL}, "\n"; #
    print "Ratio of resolved repeats in matched reads:\t", $ccounts{RR} / $ccounts{UT}, "\n";
    print "================================================================================\n";
}
