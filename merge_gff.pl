#!/usr/bin/env perl

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
my $distance   = 65535;
my $sort       = 0;
my $use_scores = 0;
my $log_scores = 0;
my $max_gap    = 1;
my $feature;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'distance|d=i' => \$distance,
    'sort|s'       => \$sort,
    'use-scores|u' => \$use_scores,
    'log-scores|l' => \$log_scores,
    'max-gap|g=i'  => \$max_gap,
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
my $previous;
print STDERR 'Merging...';
while (<$GFF_DATA>) {

    next if ($_ =~ m/^#.*$|^\s*$/);
    chomp;
    s/[\n\r]//;

    my $current = gff_read ($_);

    $current->{empty} = ($current->{score} =~ m/\d/ ? 0 : 1)
    if $use_scores;

    $current->{empty} = 1
    if !$log_scores and $current->{score} =~ m/\d/ and $current->{score} == 0;

    # if buffer has been flushed, or not initialized
    if (!defined $previous) {

        if ($current->{empty}) { # just print out window if empty
            gff_print ($current, $feature, $use_scores, $log_scores);
        }
        else { # otherwise start filling buffer, and keep track of end and start
            $previous         = $current;
            $previous->{last} = $current->{end};
        }

    }

    # adjacency rules: two windows in same chromosome, separated by at most $distance
    # and with no $max_gap adjacent empty windows are concatenated into one large window
    elsif ($previous->{seqname} =~ m/$current->{seqname}/i
           and $current->{start} - $previous->{end} <= $distance
           and ($previous->{empty} <= $max_gap)) {

        # if very first window in buffer was empty and current window is not
        # this is the actual start of the extended region
        # $previous->{first} = $current->{start} 
        # unless $current->{empty} or $previous->{first};

        # extend buffered read to current read's length
        $previous->{end} = $current->{end}
        unless $current->{end} <= $previous->{end};

        # mark this end coord as good is read is not empty
        $previous->{last} = $current->{end}
        unless $current->{empty};

        push @{$previous->{to_print}}, $current
        if $current->{empty};

        # concatenate current read's attributes to buffered read
        $previous->{attribute} .= q{; } . $current->{attribute}
        if defined $current->{attribute} and $current->{attribute} !~ m/^\.?$/;

        if ($use_scores) {
            if ($log_scores) { # convert log_2 scores to proper ones if necessary
                $previous->{score} += 2 ** $current->{score} unless $current->{empty};
            }
            else {             # and keep track of total accumulated scores
                $previous->{score} += $current->{score} unless $current->{empty};
            }
        }
        # if this is a consecutive empty window, accumulate; otherwise, reset
        $previous->{empty} = $current->{empty} ? $previous->{empty} + $current->{empty} : 0;
    }

    # print and flush buffer
    else {
        # print buffer up to the last non-empty window coordinate
        $previous->{end} = $previous->{last} if $previous->{last};
        gff_print ($previous, $feature, $use_scores, $log_scores);
        
        # print all remaining empty windows saved in buffer
        # that actually come after the last printed one
        for (@{$previous->{to_print}}) {
            gff_print ($_, $feature, $use_scores, $log_scores)
            if $_->{start} > $previous->{end}
        }

        # print the current window if it is empty and reset buffer
        if ($current->{empty}) {
            gff_print ($current, $feature, $use_scores, $log_scores);
            $previous = undef;
        }
        # or reset and start buffering windows if current is not empty
        else {
            $previous = $current;
        }
    }
}

# flush buffer if there is still something there
gff_print ($previous, $feature, $use_scores, $log_scores)
if defined $previous;

print STDERR "done\n";
## done



sub gff_print {
    my ($gff_line, $feature, $use_scores, $log_scores) = @_;

    $gff_line->{feature}   = $feature || "merged_$gff_line->{feature}";
    $gff_line->{attribute} =~ s/([;=])\s+/$1/g;
    $gff_line->{attribute} =~ s/;/; /g;

    if ($use_scores) {

        if ($log_scores) {
            $gff_line->{score}
            = sprintf ("%g", ($gff_line->{score} == 0 ? 1 : log ($gff_line->{score}) / log (2)))
        }

        $gff_line->{score} ||= q{.};

        $gff_line->{attribute} =~ s/^\.$//;
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

 merge_gff.pl [OPTION]... [FILES]...

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

 $Rev: 446 $:
 $Author: psilva $:
 $Date: 2010-12-03 14:57:39 -0800 (Fri, 03 Dec 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/merge_gff.pl $:
 $Id: merge_gff.pl 446 2010-12-03 22:57:39Z psilva $:

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
