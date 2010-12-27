#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use File::Temp qw/tempdir/;
use Cwd;

use Pod::Usage;
use Getopt::Long;

my $tmpdir = tempdir(CLEANUP => 1);
#my $tmpdir = tempdir();

my $help;
my $verbose;
my $name = 'dzlab-tools';
my $repo = 'git://dzlab.pmb.berkeley.edu/dzlab-tools-svn';
my $todos = 1;
my $tag;

my $result = GetOptions (
    "verbose|v"  => \$verbose,
    "help|h"     => \$help,
    "name|n=s"  => \$name,
    "repo|r=s"   => \$repo,
    "tag|t=s"    => \$tag,
    "todos|d"    => \$todos,
);
pod2usage(-verbose => 1) if 
(!$result || $help || !$tag || (! -d $tmpdir));  

my $dir = getcwd;
my $build = "$tmpdir/$name-$tag";
my $zipname = "$dir/$name-$tag.zip" ;

if ( ! -x '/usr/bin/todos'){
    die "please install todos utility from tofrodos";
}

if (system('git', 'clone', $repo, $build)){
    die "FATAL: can't clone";
}

chdir $build or die "FATAL: can't chdir to $build";

if (system('git', 'checkout', $tag)){
    die "FATAL: can't checkout $tag";
}
if (system('git', 'submodule', 'init')){
    die "FATAL: can't init";
}
if (system('git', 'submodule', 'update')){
    die "FATAL: can't submodule update";
}
if (system(q/find2perl -name '.git*' -prune -exec rm -rf {} \; | perl/)){
    die "FATAL: can't remove .git directories?";
}
if (system(q(find2perl -type f -exec /usr/bin/todos {} \; | perl))){
    die "FATAL: can't convert to dos line endings ?";
}
chdir $tmpdir or die "FATAL: can't chdir to $tmpdir";

if (-e $zipname && ( unlink $zipname != 1 )){
    die "FATAL: $zipname already exists but couldn't be removed first...";
}

if (system('zip', '-r', $zipname, "$name-$tag/")){
    die "FATAL: can't zip up";
}
chdir $dir; # need to get out of tmpdir for it to be removed properly

=head1 NAME

build-release.pl - build a release from a lab into a zip distribution

=head1 SYNOPSIS

build-release.pl [OPTION]... 

=head1 DESCRIPTION

Long description

=head1 OPTIONS

    --verbose     print increasingly verbose error messages
    --help        print this information

=cut

