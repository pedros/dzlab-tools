#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

# Check required command line parameters
unless (@ARGV > 0) {
    pod2usage ( -verbose => 1 );
}

# Grabs and parses command line options
my $result = GetOptions (
    'verbose|v' => sub { use diagnostics; },
    'quiet|q'   => sub { no warnings; },
    'help|h'    => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'  => sub { pod2usage ( -verbose => 2 ); }
);

my %buffer = ();
while (1) {
    my $previous = <>;
    last if !defined $previous;
    next if $previous =~ m/#|chr[cm]/i;

    my $current = <>;
    last if !defined $current;
    next if $current =~ m/#|chr[cm]/i;

    my $next = <>;
    last if !defined $next;
    next if $next =~ m/#|chr[cm]/i;

     my ($previous_c, $current_c, $next_c)
     = (0, 0, 0);

    ($previous_c) = $previous =~ m/c=(\d+)?;/;
    ($current_c)  = $current  =~ m/c=(\d+)?;/;
    ($next_c)     = $next     =~ m/c=(\d+)?;/;

    if ((exists $buffer{nnc} and $buffer{nnc}->[1] > 1) and
        ((exists $buffer{nc} and $buffer{nc} > 1) or $previous_c > 1)) {
        print $buffer{nnc}->[0];
        delete $buffer{nnc};
    }

    if ($current_c <= 1 and $previous_c > 1 and
        exists $buffer{nc} and $buffer{nc} > 1) {
        print $previous;
        delete $buffer{nc};
    }

    if ($current_c > 1) {
        print $previous if $previous_c > 1;
        print $current if $previous_c > 1 or $next_c > 1;
        print $next if $next_c > 1;
    }
    elsif ($next_c > 1) {
        my $nnext = <>;
        last if !defined $nnext;
        next if $nnext =~ m/#|chr[cm]/i;
        my ($nnext_c) = $nnext =~ m/c=(\d+)?;/;
        print $next if $nnext_c > 1;
        $buffer{nnc} = [$nnext, $nnext_c];
    }

    $buffer{pc} = $previous_c;
    $buffer{cc} = $current_c;
    $buffer{nc} = $next_c;
}


__END__


=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

=head1 REVISION

=head1 AUTHOR

 Pedro Silva <psilva@nature.berkeley.edu/>
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
