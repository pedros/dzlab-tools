#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use FindBin qw($Bin);
use File::Spec;

my $gitname = $^O eq 'MSWin32' ? 'git.exe' : 'git';

#if (! grep { -e File::Spec->catfile($_,$gitname) } File::Spec->path){
    #die "Can't find Git source code management software in PATH";
#}

chdir $Bin;

say 'u) Update DZLab-Tools to latest version';
say 'c) Checkout older version via graphical interface';
say 'q) Quit';
print 'enter u, c, or q: ';

while (my $response = <>){
    chomp $response;
    if ($response eq 'u'){
        system(qw/git fetch/);
        system(qw{git checkout remotes/origin/HEAD});
        last;
    } 
    elsif ($response eq 'c'){
        system(qw/gitk --all/);
        last;
    }
    elsif ($response eq 'q'){
        last;
    }
    else {
        say 'u) Update DZLab-Tools to latest version';
        say 'c) Checkout older version via graphical interface';
        say 'q) Quit';
        print 'enter u, c, or q: ';
    }
}

