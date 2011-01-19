#!/usr/bin/perl -na

push @ids, $F[0] unless m[^\s*\@];

END {
    %counts = map { $_ => ++$counts{$_} } @ids;
    print map  { "$_\t$counts{$_}\n" } 
          grep { not $seen{$_}++ } @ids;
}
