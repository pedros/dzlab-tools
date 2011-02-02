#!/usr/bin/env perl
# ___UNDOCUMENTED___

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'output|o=s' => \$output,
    'verbose|v'  => sub { use diagnostics; },
    'quiet|q'    => sub { no warnings; },
    'help|h'     => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'   => sub { pod2usage ( -verbose => 2 ); }
);

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

while (<>) {
    chomp;

    my %repeat = %{ read_cross_match ($_) };

    print join ("\t",
                $repeat{sequence},
                'RepeatMasker',
                'repeat',
                $repeat{start},
                $repeat{end},
                q{.},
                $repeat{strand},
                q{.},
                "Target=$repeat{id};family=$repeat{family};score=$repeat{score}",
            ),"\n";
}


sub read_cross_match {
    my $cm_line = shift;

    my @fields = split /\s/, $cm_line;

    if (@fields == 13) {
        @fields[0..7]  = @fields[0..7];
        @fields[9..13] = @fields[8..12];
        $fields[8] = q{+};
    }
    else {
        $fields[8] = q{-};
    }

    my (
        $score,
        $div,
        $del,
        $ins,
        $sequence,
        $begin,
        $end,
        $left,
        $strand,
        $id_family,
        $begin_pos,
        $end_pos,
        $begin_left,
        $linkage
    ) = @fields;

    my ($id, $family) = split /#/, $id_family;

    die $id_family if !defined $id or !defined $family;

    return {
        score    => $score,
        sequence => $sequence,
        start    => $begin,
        end      => $end,
        strand   => $strand,
        id       => $id,
        family   => $family
    };
}


__END__


=head1 NAME

 name.pl - Short description

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

 name.pl [OPTION]... [FILE]...

 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: 249 $:
 $Author: psilva $:
 $Date: 2010-01-11 21:24:34 -0800 (Mon, 11 Jan 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/cross_match_to_gff.pl $:
 $Id: cross_match_to_gff.pl 249 2010-01-12 05:24:34Z psilva $:

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
