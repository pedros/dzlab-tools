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

my $all_stats;
my $exon_stats;
my $gene_stats;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'all-stats|a=s'  => \$all_stats,
    'exon-stats|e=s' => \$exon_stats,
    'gene-stats|g=s' => \$gene_stats,
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

my %all_stats  = %{ read_stats ($all_stats,  'all' ) };
my %exon_stats = %{ read_stats ($exon_stats, 'exon') };
my %gene_stats = %{ read_stats ($gene_stats, 'gene') };

my %intron_stats     = %{ subtract_stats (\%gene_stats, \%exon_stats) };
my %intergenic_stats = %{ subtract_stats (\%all_stats,  \%gene_stats) };

my ($all_words,        $all_scores)        = print_stats (\%all_stats);
my ($exon_words,       $exon_scores)       = print_stats (\%exon_stats);
my ($gene_words,       $gene_scores)       = print_stats (\%gene_stats);
my ($intron_words,     $intron_scores)     = print_stats (\%intron_stats);
my ($intergenic_words, $intergenic_scores) = print_stats (\%intergenic_stats);

print join ("\t",
            '#word',
            '#all',
            '#exon',
            '#gene',
            '#intron',
            '#intergenic',
            '#intron/intergenic'
        ), "\n";

for my $i (0 .. @{$all_scores} - 1) {
    print $all_words->[$i], "\t";

    print join ("\t",
                $all_scores->[$i],
                $exon_scores->[$i],
                $gene_scores->[$i],
                $intron_scores->[$i],
                $intergenic_scores->[$i],
                $intron_scores->[$i] / $intergenic_scores->[$i],
            ), "\n";
}



sub print_stats {
    my ($stats) = @_;
    my $length  =  $stats->{length};
    
    my (@words, @score) = ();

    for my $word ( sort keys %{$stats->{words}} ) {
        
        my $observed = $stats->{words}{$word}{count} / $length;
        my $expected = 1;
        map {$expected *= $_ / $length} @{$stats->{words}{$word}{independent_count}};

        push @words, $word;
        push @score, $observed / $expected;

        print STDERR join ("\t",
                           $word,
                           $stats->{words}{$word}{count},
                           $length,
                           $observed,
                           $expected,
                           $observed / $expected
                       ), "\n";
    }

    return \@words, \@score;
}

sub subtract_stats {
    my ($stats_1, $stats_2) = @_;

    my %new_stats = (
        type   => "$stats_1->{type}-$stats_2->{type}",
        length => $stats_1->{length} - $stats_2->{length},
        words  => {},
    );

    for my $word ( sort keys %{$stats_1->{words}} ) {

        my @new_independent_count = ();

        for my $i (0 .. @{$stats_1->{words}{$word}{independent_count}} - 1) {

            push @new_independent_count,
            $stats_1->{words}{$word}{independent_count}->[$i] - $stats_2->{words}{$word}{independent_count}->[$i]
        }

        $new_stats{words}{$word}->{independent_count}
        = \@new_independent_count;

        $new_stats{words}{$word}->{count}
        = $stats_1->{words}{$word}{count} - $stats_2->{words}{$word}{count};

    }
    return \%new_stats;
}

sub read_stats {
    my ($stats_file, $type) = @_;

    my %stats = (
        type => $type,
        length => undef,
        words => {},
    );

    open my $STATS, '<', $stats_file or croak "Can't open $stats_file: $!";

    while (<$STATS>) {

        next if m/^\s*#|^\s*$/;
        chomp;

        my ($word, $count, $length, $observed, $expected, $likelihood, @independent_count)
        = split /\t/;

        next if $word =~ m/[^ACGT]/;

        # assert( (length $word) == scalar @independent_count,
        #    "Wrong number of independent nucleotide scores: got " . scalar @independent_count . ", but expected " . length $word);

        # assert( $likelihood == $observed / $expected,
        #        "Unexpected likelihood ratio: got $likelihood, but expected " . $observed / $expected);
        
        # assert( $observed == $count / $length,
        #         "Unexpected observed frequency: got $observed, but expected " . $count / $length);

        $stats{words}{$word} = {
            count             => $count,
            independent_count => \@independent_count
        };

        $stats{length} = $length unless defined $stats{length};
    }
    return \%stats;
}


sub assert ($;$) {
    my($assert, $name) = @_;

    unless( $assert ) {
        require Carp;
        my $msg = 'Assert failed';
        $msg .= " - '$name'" if defined $name;
        $msg .= '!';
        Carp::croak($msg);
    }
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
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/combine-nuc-freqs.pl $:
 $Id: combine-nuc-freqs.pl 249 2010-01-12 05:24:34Z psilva $:

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
