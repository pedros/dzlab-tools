#!/usr/bin/env perl
# ___UNDOCUMENTED___

use Data::Dumper;
use strict;
use warnings;
use diagnostics;

my %manifest = ();

open my $MANIFEST, '<', $ARGV[0] or die "Can't open $ARGV[0]";
while (<$MANIFEST>) {
    chomp $_;
    next if $_ =~ m/^\s*#/ or $_ =~ m/^\s*$/;
    my ($chr, $scaff, undef, undef) = split /\t/, $_;
    $manifest{$scaff} = $chr;
}
close $MANIFEST;


my %sequences = ();
my %lengths   = ();
my $chr;

open my $GENOME, '<', $ARGV[1] or die "Can't open $ARGV[1]";
while (my $line = <$GENOME>) {

    chomp $line;

    if ($line =~ m/^>([^\s]+)?\s/) {

        my $scaff = $1;

	if (exists $manifest{$scaff}) {$chr = $manifest{$scaff}}
	else {$chr = 'unknown'}

        $lengths{$chr} = 0 unless exists $lengths{$chr};
        print STDERR $chr, "\t", $scaff, "\t", $lengths{$chr}, "\n";

    }
    else {

	$sequences{$chr} .= $line;
        $lengths{$chr}   += length $line;

    }
}
close $GENOME;

exit;
foreach my $i (sort keys %sequences) {

    next if ($i =~ m/^\s*$/);
    print ">$i\n$sequences{$i}\n";

}
