#!/usr/bin/env perl

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
my $id_column = 1;
my $output;
my $missing;

# Grabs and parses command line options
my $result = GetOptions (
    'sort|r'           => \$sort,
    'list|l=s'         => \$list,
    'min-score|s=f'    => \$min_score,
    'max-score|S=f'    => \$max_score,
    'ends-tag|t=s'     => \$ends_tag,
    'id-column|c=i'    => \$id_column,
    'output|o=s'       => \$output,
    'missing|m=s'       => \$missing,
    'verbose|v'        => sub { use diagnostics; },
    'quiet|q'          => sub { no warnings; },
    'help|h'           => sub { pod2usage ( -verbose => 1 ); },
);

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

$list = index_list ($list, $ends_tag) if $list;
my $unseen = { map { $_ => 0 } keys %$list } ;
my @sorted_list;

$id_column--;

ID:
while (<>) {
    chomp;

    my $id = (split /\t/, $_)[$id_column];
    $id =~ s/$ends_tag//;

    next ID unless
    exists $list->{$id}
    and $list->{$id}->[0] ne q{.}
    and ($min_score == 0 or $list->{$id}->[0] >= $min_score)
    and ($max_score == 0 or $list->{$id}->[0] <= $max_score);

    delete $unseen->{$id};

    if ($sort) {
        $sorted_list[$list->{$id}->[2]] = $_;
    }
    else {
        print "$_\n";
    }
}

print map { "$_\n" } grep {defined $_} @sorted_list if $sort;

if ($missing){
    use autodie;
    open my $m, '>', $missing;
    foreach my $miss (keys %$unseen) {
        print $m "$miss\n";
    }
    close $m;
}


sub index_list {
    my ($list, $ends_tag) = @_;

    my $counter = 0;
    my %list = ();
    open my $LIST, '<', $list or croak "Can't open $list: $!";
    while (<$LIST>) {

        my ($id, $freq, $alt) = (0, 0, 0);

        my @fields = split /\t/, $_;

        if (@fields > 1 and @fields < 9) {
            ($id, $freq, $alt) = @fields;
        }
        elsif (@fields > 8) {
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
 -c, --id-column   select which column in file to filter contains the locus id (default: 1)
 -m, --missing     filename to write ID from list which are not present in file
 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information

=head1 REVISION

 Version 0.0.1

 $Rev: 375 $:
 $Author: psilva $:
 $Date: 2010-07-02 16:37:53 -0700 (Fri, 02 Jul 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/filter_ends.pl $:
 $Id: filter_ends.pl 375 2010-07-02 23:37:53Z psilva $:

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
