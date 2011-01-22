#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use XML::Simple;
use XML::Twig;

my $xml_hits;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'xml-hits|x=s' => \$xml_hits,
    'output|o=s' => \$output,
    'verbose|v'  => sub { use diagnostics; },
    'quiet|q'    => sub { no warnings; },
    'help|h'     => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'   => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless $xml_hits or $xml_hits = $ARGV[0];

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}


# create object
my $xml = new XML::Simple(
    NormaliseSpace      => 2,
    ForceArray          => 0,
    KeyAttr             => {
        begin => 'content',
        end   => 'content'
    },
    NoAttr              => 0,
    ContentKey          => '-content',
);

my $t = XML::Twig->new( 
    # the twig will include just the root and selected titles 
    twig_roots => {'segments' => \&print_n_purge},
);

$t->parsefile($xml_hits);


sub print_n_purge {
    my( $t, $elt)= @_;
    my $tree = $xml->XMLin($t->sprint);

    for my $node ( $tree->{segments} ) {

        my $chromosome = $node->{query};

        next unless ref $node->{segment} eq 'ARRAY';

        for my $segment (@{$node->{segment}}) {

            next unless defined $segment;

            my $source       = $segment->{source};
            my $significance = $segment->{significance};
            my $name         = $segment->{name};
            my $score        = $segment->{score};
            my $begin        = $segment->{begin};
            my $end          = $segment->{end};
            my $number       = $segment->{segment_number};
            
            print join ("\t",
                        $chromosome,
                        'repeat_runner',
                        'repeat',
                        $begin,
                        $end,
                        $score,
                        q{.},
                        q{.},
                        "ID=$chromosome.$number; Name=$number; Note=$name; e=$significance",
                    ), "\n";
        }
    }

    $t->purge;
}


__END__


=head1 NAME

 xml_to_gff.pl - Convert RepeatRunner's XML format to GFF3

=head1 SYNOPSIS

 xml_to_gff.pl -x repeat-runner_out.xml -o repeat-runner_out.gff

=head1 DESCRIPTION

=head1 OPTIONS

 xml_to_gff.pl [OPTION]... [FILE]...

 -x. --xml         input filename in RepeatRunner's XML format
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
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/xml_to_gff.pl $:
 $Id: xml_to_gff.pl 249 2010-01-12 05:24:34Z psilva $:

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
