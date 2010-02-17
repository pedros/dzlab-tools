#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long; Getopt::Long::Configure ('bundling');
use Pod::Usage;

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my $list;
my $min_score = 0;
my $max_score = 0;
my $ends_tag  = 'ID=';
my $sort      = 0;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'sort|r'           => \$sort,
    'list|l=s'         => \$list,
    'min-score|s=f'    => \$min_score,
    'max-score|S=f'    => \$max_score,
    'ends-tag|t=s'     => \$ends_tag,
    'output|o=s'       => \$output,
    'verbose|v'        => sub { use diagnostics; },
    'quiet|q'          => sub { no warnings; },
    'help|h'           => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'         => sub { pod2usage ( -verbose => 2 ); }
);

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

$list = index_list ($list) if $list;
my @sorted_list;

ID:
while (<>) {
    my ($id, undef) = split /\t/;

    $id =~ s/$ends_tag//;

    next ID unless
    exists $list->{$id}
    and $list->{$id}->[0] ne q{.}
    and ($min_score == 0 or $list->{$id}->[0] >= $min_score)
    and ($max_score == 0 or $list->{$id}->[0] <= $max_score);

    if ($sort) {
        $sorted_list[$list->{$id}->[2]] = $_;
    }
    else {
        print $_;
    }
}

print grep {defined $_} @sorted_list if $sort;


sub index_list {
    my ($list) = @_;

    my $counter = 0;
    my %list = ();
    open my $LIST, '<', $list or croak "Can't open $list for reading";
    while (<$LIST>) {

        my ($id, $freq, $alt) = (0, 0, 0);

        my @fields = split /\t/, $_;

        if (@fields < 9) {
            ($id, $freq, $alt) = @fields;
        }
        elsif (@fields == 9) {
            ($id, $freq, $alt) = @fields[8, 5, 0];
            $id =~ s/^.*$ends_tag([^;]+).*$/$1/;
        }
        else {
            ($id) = $fields[0];
        }

        $id =~ s/[\r\n]//g;
        
        $list{$id} = [$freq, $alt, $counter++];
    }
    close $LIST or carp "Can't close $list after reading";

    return \%list;
}


__END__


=head1 NAME

 filter_ends.pl - Filter an ends analysis output file by a list of gene IDs

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

 filter_ends.pl [OPTION]... [FILE]...

 -l, --list        filename with one or two fields: ID (required) and score (optional)
 -s, --min-score   minimum score to filter by (0 by default)
 -S, --max-score   maximum score to filter by (0 by default)
 -r, --sort        sort input file by order in list
 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev$:
 $Author$:
 $Date$:
 $HeadURL$:
 $Id$:

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
