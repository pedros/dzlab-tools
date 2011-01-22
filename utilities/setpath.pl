#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use File::Find;
use File::Spec::Functions;

my $dest = '/usr/local/bin/';

find(sub {
        # $File::Find::dir  = /some/path/
        # $_                = foo.ext
        # $File::Find::name = /some/path/foo.ext
        if (-x $File::Find::name){
            symlink($File::Find::name,catfile($dest,$_));
        }
    }, qw(
    /cygdrive/c/strawberry/c/bin
    /cygdrive/c/strawberry/perl/site/bin
    /cygdrive/c/strawberry/perl/bin
    ));
