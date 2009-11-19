#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

my $DATA_HANDLE = 'ARGV';
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'output|o=s'  => \$output,
    'verbose|v'   => sub { use diagnostics; },
    'quiet|q'     => sub { no warnings; },
    'help|h'      => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'    => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my %contexts = (
    x => 'CG',  X => 'CG',
    y => 'CHG', Y => 'CHG',
    z => 'CHH', Z => 'CHH'
);
my %c_map;

while (<$DATA_HANDLE>) {
    chomp;
    my $seeker = parse_seeker ($_);

    update_counts ($seeker, \%c_map, \%contexts);
    
}

sum_counts (\%c_map);


sub parse_seeker {
    my ($seeker_read) = @_;
    
    my @seeker_fields = split /[\t\s]/, $seeker_read;

    return undef unless @seeker_fields == 8;

    my %seeker = (
        read_id       => $seeker_fields[0],
        mismatches    => $seeker_fields[1],
        read_strand   => ($seeker_fields[2] =~ m{\+FW|-RC} ? q{+} : q{-}),
        ref_strand    => substr ($seeker_fields[3], 2, 1),
        chromosome    => substr ($seeker_fields[3], 0, 2),
        coordinate    => substr ($seeker_fields[3], 3, 10),
        ref_sequence  => (split /_/, $seeker_fields[4])[1],
        read_sequence => $seeker_fields[5],
        summary       => $seeker_fields[6],
        index         => $seeker_fields[7],
    );

    return \%seeker;
}


sub update_counts {
    my ($seeker, $c_map, $contexts) = @_;

    for my $idx (0 .. (length $seeker->{summary}) - 1) {
        my $base = chr vec ($seeker->{summary}, $idx, 8);
        next if $base eq q{-};

        my $chr = $seeker->{chromosome};           # chromosome
        my $con = $contexts->{$base};              # context
        my $crd = $seeker->{coordinate} + $idx;    # coordinate
        my $met = $base =~ m/[xyz]/ ? q{t} : q{c}; # c or t
        
        $c_map->{$chr}{$con}{$crd}{$met}++;        # increase c or t count
        $c_map->{$chr}{$con}{$crd}{s}              # save strand
            = $seeker->{ref_strand};
    }
}

sub sum_counts {
    my ($c_map) = @_;

    for my $chr (sort keys %$c_map) {
        for my $con (sort keys %{$c_map->{$chr}}) {
            for my $crd (sort {$a <=> $b} keys %{$c_map->{$chr}{$con}}) {

                my $c = $c_map->{$chr}{$con}{$crd}{c} || 0;
                my $t = $c_map->{$chr}{$con}{$crd}{t} || 0;
                my $s = $c_map->{$chr}{$con}{$crd}{s} || q{+};

                my $methylation
                    = ($c + $t > 0 ? sprintf ("%g", $c / ($c + $t)) : 'NaN');
                
                print join ("\t",
                            $chr,
                            'bsseeker',
                            $con,
                            $crd,
                            $crd,
                            $methylation,
                            $s,
                            q{.},
                            "c=$c; t=$t",
                        ), "\n";
            }
        }
    }
}



__END__


=head1 NAME

 seeker2gff.pl - Convert BS Seeker's output format to GFF version 3

=head1 SYNOPSIS

 seeker2gff.pl data.seeker -o data.gff

=head1 DESCRIPTION

=head1 OPTIONS

 seeker2gff.pl [OPTION]... [FILE]...

 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: $:
 $Author: $:
 $Date: $:
 $HeadURL: $:
 $Id: $:

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
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
