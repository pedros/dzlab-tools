#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my $GFF_DATA   = 'ARGV';
my $distance   = 50;
my $sort       = 0;
my $use_scores = 0;
my $log_scores = 0;
my $feature;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'distance|d=i' => \$distance,
    'sort|s'       => \$sort,
    'use-scores|u' => \$use_scores,
    'log-scores|l' => \$log_scores,
    'feature|f=s'  => \$feature,
    'output|o=s'   => \$output,
    'verbose|v'    => sub { use diagnostics; },
    'quiet|q'      => sub { no warnings; },
    'help|h'       => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'     => sub { pod2usage ( -verbose => 2 ); }
);

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

if ($sort) {
    print STDERR 'Sorting...';
    my $sorted_filename = sort_repeats (@ARGV);
    $GFF_DATA = undef;
    open $GFF_DATA, '<', $sorted_filename;
    print STDERR "done\n";
}

# one-step buffer
my $previous = undef;

print STDERR 'Merging...';
while (<$GFF_DATA>) {

    next if ($_ =~ m/^#.*$|^\s*$/);
    chomp;
    s/[\n\r]//;

    my $current = gff_read ($_);

    # check empty windows
    $current->{empty} = ($current->{score} =~ m/\d/ ? 0 : 1)
    if $use_scores;

    $current->{empty} = 1
    if !$log_scores and $current->{score} == 0;

    # if buffer has been flushed, or not initialized
    if (!defined $previous) {
        $previous = $current unless $current->{empty};
        $previous->{last} = $current->{end} unless $current->{empty};
    }

    # adjacency rules: two windows in same chromosome, separated by at most $distance
    # and with no two adjacent empty windows, are concatenated into one large window
    elsif ($previous->{seqname} =~ m/$current->{seqname}/i
           and $current->{start} - $previous->{end} <= $distance
           and ($use_scores == 0 || $previous->{empty} == 0 or $current->{empty} == 0)) {

        $previous->{start} = $current->{start}
        if $previous->{empty};

        # extend buffered read to current read's length
        $previous->{end} = $current->{end}
        unless $current->{end} <= $previous->{end};

        # concatenate current read's attributes to buffered read
        $previous->{attribute} .= q{; } . $current->{attribute}
        if defined $current->{attribute} and $current->{attribute} ne q{.}
        and $current->{attribute} ne q{};

        if ($use_scores) {
            if ($log_scores) {
                $previous->{total} += 2 ** $current->{score} unless $current->{empty};
            }
            else {
                $previous->{total} += $current->{score} unless $current->{empty};
            }
            $previous->{score}  = $current->{score};
            $previous->{empty}  = $current->{empty};
            $previous->{last}   = $current->{end} unless $current->{empty};
        }
    }

    # print and flush buffer
    else {
        $previous->{end} = $previous->{last} if $previous->{last};
        gff_print ($previous, $feature, $use_scores, $log_scores);
        $previous = undef;
        $previous = $current unless $previous->{empty};
    }
}

gff_print ($previous, $feature) if defined $previous;

print STDERR "done\n";
## done



sub gff_print {
    my ($gff_line, $feature, $use_scores, $log_scores) = @_;


    $gff_line->{feature}   = $feature || "merged_$gff_line->{feature}";
    $gff_line->{attribute} =~ s/([;=])\s+/$1/g;
    $gff_line->{attribute} =~ s/;/; /g;

    if ($use_scores) {
        $gff_line->{total} = 0 unless defined $gff_line->{total};

        if ($log_scores) {
            $gff_line->{score}
            = sprintf ("%g", ($gff_line->{total} == 0 ? 1 : log ($gff_line->{total}) / log (2)))
        }
        else {$gff_line->{score} = $gff_line->{total}}

        $gff_line->{score} = q{.} unless $gff_line->{total};

        $gff_line->{attribute} = s/^\.$//;
        $gff_line->{attribute} .= 'length=' . ($gff_line->{end} - $gff_line->{start} + 1);
    }

    print join ("\t",
                $gff_line->{seqname},
                $gff_line->{source},
                $gff_line->{feature},
                $gff_line->{start},
                $gff_line->{end},
                $gff_line->{score},
                $gff_line->{strand},
                $gff_line->{frame},
                $gff_line->{attribute},
            ), "\n";
}        


sub gff_read {
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute)
    = split /\t/, shift;

    return {
	'seqname'  => $seqname,
	'source'   => $source,
	'feature'  => $feature,
	'start'    => $start,
	'end'      => $end,
	'score'    => $score,
	'strand'   => $strand,
	'frame'    => $frame,
	'attribute'=> $attribute
    };
}


sub sort_repeats {
    use File::Temp qw/tempfile/;

    my @repeats_files           = @_;
    my @repeats                 = ();
    my ($tmp_fh, $tmp_filename) = tempfile();

    BATCH:
    for my $repeats_file (@repeats_files) {

        open my $GFF, '<', $repeats_file 
        or croak "Can't read $repeats_file: $!";

        while (<$GFF>) {
            next if ($_ =~ m/^#.*$|^\s*$/);
            chomp;
            push @repeats, $_;
        }

        close $GFF
        or carp "Can't close $repeats_file: $!";
    }

    map  { print $tmp_fh $_, "\n" }
    sort {
        (split /\t/, $a)[0] cmp (split /\t/, $b)[0] or
        (split /\t/, $a)[3] <=> (split /\t/, $b)[3]
    }
    @repeats;

    return $tmp_filename;
}


__END__


=head1 NAME

 merge_gff.pl - Merge multiple GFF files based on coordinates

=head1 SYNOPSIS

 merge_gff.pl --distance 50 --sort file1.gff file2.gff ... filen.gff

=head1 DESCRIPTION

=head1 OPTIONS

 merge_gff.pl [OPTION]... [FILE]...

 -d, --distance    maximum distance between adjacent features to be merged
 -s, --sort        sort input files by sequence name and starting coordinate
 -u, --use-scores  look for empty windows and use two adjacent ones to break merging
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
