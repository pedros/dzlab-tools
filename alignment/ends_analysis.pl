#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(min max);

my $DATA_HANDLE    = 'ARGV';
my $gff_annotation = q{};
my $bin_width      = 100;
my $distance       = 5000;
my $stop_flag      = 2;
my $stop_distance  = 1500;
my $three_prime    = 0;
my $five_prime     = 0;
my $attribute_id   = 'ID';
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'gff-annotation|g=s' => \$gff_annotation,
    'bin-width|b=i'      => \$bin_width,
    'distance|d=i'       => \$distance,
    'stop-flag|s=i'      => \$stop_flag,
    'stop-distance|k=i'  => \$stop_distance,
    'three-prime|3'      => \$three_prime,
    'five-prime|5'       => \$five_prime,
    'extract-id|x=s'     => \$attribute_id,
    'output|o=s'         => \$output,
    'verbose|v'          => sub { use diagnostics; },
    'quiet|q'            => sub { no warnings; },
    'help|h'             => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'           => sub { pod2usage ( -verbose => 2 ); }
);

my $stop_flag_dispatch = {
    0 => \&stop_flag_0,
    1 => \&stop_flag_1,
    2 => \&stop_flag_2,
    3 => \&stop_flag_3,
    4 => \&stop_flag_4,
    5 => \&stop_flag_5,
    6 => \&stop_flag_6,
};

# Check required command line parameters
pod2usage ( -verbose => 1 )
    unless @ARGV and $result
    and $gff_annotation
    and ($three_prime xor $five_prime)
    and exists $stop_flag_dispatch->{$stop_flag};

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't read $output: $!";
    select $USER_OUT;
}

print STDERR 'Indexing and offsetting annotation file...';
my $annotation
    = offset_gff_annotation (index_gff_annotation ($gff_annotation, $attribute_id),
                             $stop_flag_dispatch->{$stop_flag},
                             {
                                 -three_prime   => $three_prime,
                                 -five_prime    => $five_prime,
                                 -distance      => $distance,
                                 -stop_flag     => $stop_flag,
                                 -stop_distance => $stop_distance
                             });
print STDERR "done\n";

my $gff_it = make_gff_iterator ($ARGV[0], \&gff_read);

my $num_bins = 2 * int ($distance / $bin_width);
my @totals_bins;

my %sorted_annotations;

print STDERR 'Assigning coordinates to bins...';
while (my $gff_line = $gff_it->()) {

    next unless ref $gff_line eq 'HASH';
    
    print STDERR $gff_line->{start}, "\b" x length $gff_line->{start};
    my $brs = binary_range_search (
        [ $gff_line->{start}, $gff_line->{end} ],
        $sorted_annotations{$gff_line->{seqname}}
            ||= [ map { $annotation->{$gff_line->{seqname}}{$_} }
                      sort { $a <=> $b }
                          keys %{$annotation->{$gff_line->{seqname}}}]
    );
    
    next unless $brs;

    my $reverse
        = ($five_prime and $gff_line->{strand} eq q{-})
            or ($three_prime and $gff_line->{strand} eq q{+}) || 0;

    my $index
        = abs int (($brs->[$reverse] - ( int ($gff_line->{end} - $gff_line->{start}) / 2 + $gff_line->{start})) / $bin_width);

    $index = $reverse ? $num_bins - 1 - $index : $index;
    
    push @{$totals_bins[$index]}, $gff_line->{score};
}
print STDERR "done\n";


print STDERR 'Computing scores...';
for my $index (0 .. $num_bins - 1) {

    my $total_score;
    unless (defined $totals_bins[$index]) {
        $totals_bins[$index] = 'NaN';
    }
    else {
        for my $score (0 .. @{$totals_bins[$index]} - 1) {
            $total_score += $totals_bins[$index]->[$score];
        }
        $totals_bins[$index] = $total_score / @{$totals_bins[$index]};
    }
    print $index * $num_bins - $distance, "\t", $totals_bins[$index], "\n";
}
print STDERR "done\n";

 
sub binary_range_search {
    my ($range, $ranges) = @_;
    my ($low, $high)     = (0, @{$ranges} - 1);
    while ($low <= $high) {
        my $try = int (($low + $high) / 2);
        $low  = $try + 1, next if $ranges->[$try][1] < $range->[0];
        $high = $try - 1, next if $ranges->[$try][0] > $range->[1];
        return $ranges->[$try];
    }
    return;
}

sub make_gff_iterator {
    my ($gff_file_name, $gff_parser) = @_;
    
    open my $GFF_FILE, '<', $gff_file_name or croak "Can't read $gff_file_name: $!";
    
    return sub { $gff_parser->(scalar <$GFF_FILE>) };
}

sub gff_read {
    return [] if $_[0] =~ m/^\s*#+/;

    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute)
    = split /\t/, shift || return;

    $attribute =~ s/[\r\n]//g;

    return {
	'seqname'   => lc $seqname,
	'source'    => $source,
	'feature'   => $feature,
	'start'     => $start,
	'end'       => $end,
	'score'     => $score,
	'strand'    => $strand,
	'frame'     => $frame,
	'attribute' => $attribute
    };
}

sub index_gff_annotation {
    my ($gff_file_name, $attribute_id) = @_;

    my $gff_annotation = {};

    my $gff_it = make_gff_iterator ($gff_file_name, \&gff_read);
    
    while (my $gff_line = $gff_it->()) {

        $gff_line->{attribute} =~ s/.*$attribute_id[=\s]?([^;\s]+).*/$1/
            if $attribute_id;
        
        $gff_annotation->{$gff_line->{seqname}}{$gff_line->{start}}
            = [$gff_line->{start}, $gff_line->{end}, $gff_line->{strand}, $gff_line->{attribute}];
    }

    return $gff_annotation;
}

sub offset_gff_annotation {
    my ($gff_annotation, $flag_parser, $parameters) = @_;

    return unless $gff_annotation and $flag_parser;
    
    for my $seqid (sort keys %{$gff_annotation}) {

        my @memory = ([0, 0, q{+}, q{}]);

        for my $start (sort {$a <=> $b} keys %{$gff_annotation->{$seqid}}) {

            push @memory, $gff_annotation->{$seqid}{$start} if @memory < 3;

            if (@memory == 3) {
                my ($flag_start, $flag_end)
                    = $flag_parser->($parameters, @memory);

                delete $gff_annotation->{$seqid}{$start};
                $gff_annotation->{$seqid}{$flag_start}
                    = [$flag_start, $flag_end, $memory[1]->[2], $memory[1]->[3]];
                shift @memory;
            }
        }
    }
    return $gff_annotation;
}


sub stop_flag_2 {
    my ($parameters, $previous, $current, $next) = @_;

    my ($prev_end, $next_start) = ($previous->[1], $next->[0]);
    my ($minimum, $maximum)     = min_max_distance ($current, $parameters);

    my $flag_start = max ($prev_end, $minimum);
    my $flag_end   = min ($next_start, $maximum);
    
    return ($flag_start, $flag_end);
}


sub min_max_distance {
    my ($current, $parameters) = @_;

    my ($minimum, $maximum);

    if ($parameters->{-five_prime} and $current->[2] eq q{-}
            or $parameters->{-three_prime} and $current->[2] eq q{+}) {

        $minimum = $current->[1] - $parameters->{-distance};
        $maximum = $current->[1] + $parameters->{-distance};

    }
    else {
        $minimum = $current->[0] - $parameters->{-distance};
        $maximum = $current->[0] + $parameters->{-distance};
    }

    $minimum = ($minimum < 0 ? 0 : $minimum);

    return ($minimum, $maximum);
}


sub stop_flag_0 {
    croak "Not yet implemented";
}

sub stop_flag_1 {
    croak "Not yet implemented";
}


sub stop_flag_3 {
    croak "Not yet implemented";
}

sub stop_flag_4 {
    croak "Not yet implemented";
}

sub stop_flag_5 {
    croak "Not yet implemented";
}

sub stop_flag_6 {
    croak "Not yet implemented";
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
