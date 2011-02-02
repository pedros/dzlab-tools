#!/usr/bin/env perl
# ___UNDOCUMENTED___

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/sum/;

use FindBin;
use lib "$FindBin::Bin/DZLab-Tools/lib";
use DZLab::Tools::Fasta;

my $GFF_DATA  = 'ARGV';
my $min_score = 0.10;
my $min_total = 10;
my $min_gap   = 1;
my $normalize = q{};
my $sort      = 0;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'min-score|s=f' => \$min_score,
    'min-total|t=i' => \$min_total,
    'min-gap|g=i'   => \$min_gap,
    'normalize|n=s' => \$normalize,
    'sort|r'        => \$sort,
    'output|o=s'  => \$output,
    'verbose|v'   => sub { use diagnostics; },
    'quiet|q'     => sub { no warnings; },
    'help|h'      => sub { pod2usage ( -verbose => 1 ); },
    'manual|p'    => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result and $min_score > 0;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

if ($sort) {
    my $sorted_filename = sort_gff (@ARGV);
    $GFF_DATA = undef;
    open $GFF_DATA, '<', $sorted_filename;
}

my $reference = slurp_fasta ($normalize, {-l => 1});

my %buffer = ();
my %counts = ();

while (<$GFF_DATA>) {
    next if m/^#.*$|^\s*$/;
    chomp;
    s/[\n\r]//;

    my $current = gff_read ($_);
    my $seqid   = $current->{seqname};
    my $score   = $current->{score};
    my $total   = sum map { (split /=/, $_)[1] } split /;/, $current->{attribute};

    # initialize buffer
    if (!exists $buffer{$seqid} or !@{$buffer{$seqid}}) {
        push @{$buffer{$seqid}}, [$score, $total, $current->{end}]
        if $score ne q{.} and $score > $min_score and $total > $min_total;
    }
    # extend buffer if possible
    elsif ($current->{start} - $buffer{$seqid}->[-1]->[2] <= $min_gap
           and $score ne q{.} and $score > $min_score and $total > $min_total) {
        push @{$buffer{$seqid}}, [$score, $total, $current->{end}];
    }
    # flush buffer
    else {
        flush_buffer (\%buffer, \%counts, $seqid);
    }
}

my %full_counts = ();

# post-process histogram counts, one chromosome at a time
for my $seqid (keys %counts) {

    # flush whatever's left over, if anything
    flush_buffer (\%buffer, \%counts, $seqid)
    if @{$buffer{$seqid}};
 
    # slide through each bin size in ascending sorted order
    for my $size (sort { $a <=> $b } keys %{$counts{$seqid}}) {

        # find number of hits per bin, mean score and deviations
        my $number = scalar @{$counts{$seqid}{$size}{mean}};
        # my $mean   = (sum @{$counts{$seqid}{$size}{mean}}) / $number;
        # my $var    = (sum map { ($_ - $mean) ** 2 }
        #               @{$counts{$seqid}{$size}{mean}}) / ($number - 1 || 1);
        # my $std    = $var ** 1/2;
    
        # normalize number of hits with bulk size of genome
        if ($normalize and $reference and exists $reference->{$seqid}) {
            $number /= (length $reference->{$seqid});
            $number  = sprintf ("%g", $number * 1000000);
        }

        # print join ("\t",
        #             $size,
        #             $number,
        #             sprintf ("%g", $mean),
        #             sprintf ("%g", $var),
        #             sprintf ("%g", $std),
        #         ), "\n";

        $full_counts{$size} += $number;
    }
}

#print "#$seqid, counts normalized per 1000bp\n";
#print "#Window size\tWindow number\tMean score\tVariance\tStd\n";

for my $size (sort {$a <=> $b } keys %full_counts) {
    print $size, "\t", $full_counts{$size}, "\n";
}


sub flush_buffer {
    # alters $buffer and $counts in place
    my ($buffer, $counts, $seqid) = @_;

    my $size = scalar @{$buffer->{$seqid}};
    my $mean = (sum map {  $_->[0]               } @{$buffer->{$seqid}}) / $size || 1;
    # my $var  = (sum map { ($_->[0] - $mean) ** 2 } @{$buffer->{$seqid}}) / ($size - 1 || 1);

    push @{$counts->{$seqid}{$size}{mean}}, $mean;
    # push @{$counts->{$seqid}{$size}{var }}, $var;
    # push @{$counts->{$seqid}{$size}{std }}, $var ** 1/2;
    
    @{$buffer->{$seqid}} = ();
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


sub sort_gff {
    use File::Temp qw/tempfile/;

    my @gff_files               = @_;
    my @gff                     = ();
    my ($tmp_fh, $tmp_filename) = tempfile();

    BATCH:
    for my $gff_file (@gff_files) {

        open my $GFF, '<', $gff_file 
        or croak "Can't read $gff_file: $!";

        while (<$GFF>) {
            next if ($_ =~ m/^#.*$|^\s*$/);
            chomp;
            push @gff, $_;
        }

        close $GFF
        or carp "Can't close $gff_file: $!";
    }

    map  { print $tmp_fh $_, "\n" }
    sort {
        (split /\t/, $a)[0] cmp (split /\t/, $b)[0] or
        (split /\t/, $a)[3] <=> (split /\t/, $b)[3]
    }
    @gff;

    return $tmp_filename;
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
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/count_adjacencies.pl $:
 $Id: count_adjacencies.pl 249 2010-01-12 05:24:34Z psilva $:

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
