#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Spec;
use File::Path;
use File::Basename;

my $name;
my $genome;
my $contexts;
my $width = 50;
my $step = 50;
my $remove;
my $work_dir = File::Spec->curdir();
my $mode = 'overwrite';
my %mode_table = (      # might do this later
    'overwrite' => "", 
    'append' => "");

my $result = GetOptions(
    'name|n=s'     => \$name,
    'genome|g=s'   => \$genome,
    'contexts|c=s' => \$contexts,
    'width|w=i'    => \$width,
    'step|s=i'     => \$step,
    'work_dir|d=s' => \$work_dir,
    'mode|m=s'     => \$mode
    );
    


# run_cmd ("cd $NAME");  


croak "Bad directory"
  unless -e "$work_dir/post-processing/single-c"
      and -e "$work_dir/post-processing/windows";


if ($mode == 'overwrite') {
    run_cmd ("rm -r post-processing/single-c");

}



print "Concatenating single c files";


















sub run_cmd {
    my ($cmd) = @_;
    warn "-- CMD: $cmd\n";
    eval {system ("$cmd")};
    croak "** failed to run command '$cmd': $@" if $@;
}
















=head1 NAME

 cat_and window.pl - concatenate, merge, and window single c files

=head1 SYNOPSIS


=head1 DESCRIPTION

=head1 OPTIONS

 cat_and_window-pl [OPTION]... [FILE]...

 -m,   --mode
 -n,   --name
 -g,   --genome
 -c,   --contexts
 -w,   --width
 -s,   --step
 -d,   --working-directory

=head1 REVISION

 

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
 Jon Graff <jonthegraff@berkeley.edu>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
