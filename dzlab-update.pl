#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use IO::Prompt;
use FindBin qw($Bin);
use File::Spec;

my $gitname = $^O eq 'MSWin32' ? 'git.exe' : 'git';

if (! grep { -e File::Spec->catfile($_,$gitname) } File::Spec->path){
    die "Can't find Git source code management software in PATH";
}

chdir $Bin;

my $response1 = prompt -menu=>[
'Update DZLab-Tools to latest version',
'Checkout older version via graphical interface',
'Quit'
];

if ($response1 =~ /^Update/){
    system(qw/git fetch/);
    system(qw{git checkout remotes/origin/HEAD});
} 
elsif ($response1 =~ /^Checkout/){
    system(qw/gitk --all/);
}


