#!/usr/bin/env perl
# ___UNDOCUMENTED___

my %hash = ();
while (<>) {
    my @a=split (/\t/, $_);
    $hash{"$a[1]"}++;
}

my $total = 0;
foreach my $key (keys %hash) {
    $total += $hash{$key};
}

foreach my $key (sort {$a <=> $b} keys %hash) {
    print $key, "\t", $hash{$key}, "\t", sprintf("%.3f", $hash{$key}/$total*100), "%\n";
}

print "Total\t$total\t100%\n";
