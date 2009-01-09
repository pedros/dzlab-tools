#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;disable diagnostics;
use Data::Dumper;
use Getopt::Long;
use Carp;
use Pod::Usage;

# Globals, passed as command line options
my $gff_file = q{-};
my $width = 0;
my $step = 0;
my $no_sort = 0;
my @batch = ();
my $output = q{-};
my $verbose = 0;
my $quiet = 0;

my @argv = @ARGV;

# Grabs and parses command line options
my $result = GetOptions (
    'gff-file|f=s' => \$gff_file,
    'width|w=i' => \$width,
    'step|s=i' => \$step,
    'no-sort|n' => \$no_sort,
    'output|o:s' => \$output,
    'batch|b=s{,}' => \@batch,
    'verbose|v' => sub {enable diagnostics;use warnings;},
    'quiet|q' => sub {disable diagnostics;no warnings;},
    'help|usage|h' => sub {pod2usage(-verbose => 1);},
    'manual|man|m' => sub {pod2usage(-verbose => 2);}
    );

# Check required command line parameters
if ( ! $width ) { pod2usage(-verbose => 1) }

BATCH:
while (@batch or $gff_file) {
    if (@batch) {
        print STDERR "Running in batch mode: all input files will have \'.avg\' appended to form the output name.\n";
        $gff_file = shift @batch;
        $output = $gff_file . "_w${width}_s${step}.avg";
    }

    # redirects STDOUT to file if specified by user
    if ($output ne '-') {
        open STDOUT, '>', "$output" or croak "Can't redirect STDOUT to file: $output";
    }

    # opens gff file or STDIN
    my $GFF;
    if ($gff_file ne '-') {open $GFF, '<', $gff_file or croak "Can't read file: $gff_file"}
    else {$GFF = 'STDIN'}

    # reads in data
    my @data;
    while (<$GFF>) {
        chomp;
        next if ($_ =~ m/^#.*$|^\s*$/);
        push @data, $_;
    }
    if ($GFF ne 'STDIN') {close $GFF}

    # prints out header fields that contain gff v3 header, generating program, time, and field names
    gff_print_header ($0, @argv);

    # sort databy sequence name, feature, and starting coordinates
    unless ($no_sort) {
        @data =  gff_sort ( \@data );
    }

    # gets indices for chromosome changes (tracks changes on 0th field - seqname)
    my @locsindex = gff_find_array ( 0, \@data );

    # for each chromosome
    for (my $i = 0; $i < @locsindex; $i++) {

        # get an array slice with only that chromosome
        my @subdata;
        if ($i + 1 <@locsindex) {
            @subdata = @data[$locsindex[$i]..$locsindex[$i+1]];
        } else {
            @subdata = @data[$locsindex[$i]..$#data]; # prevents array pointer out of bounds
        }

        # gets indices for context changes (tracks changes on 2th field - feature)
        my @featureindex = gff_find_array ( 2, \@subdata );

        # for each feature/context
        for (my $l = 0; $l < @featureindex; $l++) {

            # get an array slice containing just that feature/chromosome and feed it to the sliding window routine
            if ($l + 1 < @featureindex) {
                gff_sliding_window ( $width, $step, @subdata[$featureindex[$l]..$featureindex[$l + 1] - 1]);
            } else {
                gff_sliding_window ( $width, $step, @subdata[$featureindex[$l]..$#subdata]);
            }
        }
    }
    $gff_file = 0;
}

close STDOUT;
exit 0;

sub gff_sliding_window {
    my ($width, $step, @data) = @_;

    my %last_rec = %{gff_read ($data[-1])};

    print STDERR
	"windowing data on sequence $last_rec{'seqname'} and context $last_rec{'feature'} with $width bp width and $step bp step...\n";

    my $lastcoord = $last_rec{'end'};
    my $seqname = $last_rec{'seqname'};
    my $context = $last_rec{'feature'};
    my $lastrecord = 0;

    for (my $i = 1; $i < $lastcoord; $i++) {

	my ($c_count, $t_count, $score) = (0, 0, -0.1);
        my ($overlap_c, $overlap_t) = (0, 0);

	my @range = @{ gff_filter_by_coord ($i, $i + $width, $lastrecord, \@data) };

	$lastrecord = shift @range;

	foreach my $k (@range) {
	    my %current_rec = %{gff_read ($k)};
	    my ($c_tmp, $t_tmp) = split(/;/, $current_rec{'attribute'});
	    ($c_tmp) = $c_tmp =~ m/(\d+)/;
	    ($t_tmp) = $t_tmp =~ m/(\d+)/;
            $c_count += $c_tmp;
            $t_count += $t_tmp;
	    $seqname = $current_rec{'seqname'};
	    $context = $current_rec{'feature'};

            if ($current_rec{'start'} > $i + $step) {
                $overlap_c += $c_tmp;
                $overlap_t += $t_tmp;
            }
	}

	if ($c_count + $t_count != 0) {
	    if ($c_count != 0) {$score = $c_count / ($c_count + $t_count)}
	}

	my $attribute = "c=$c_count;t=$t_count";
        if ($step) {$attribute .= ";over_c=$overlap_c;over_t=$overlap_t"}

	if (scalar(@range) == 0) {
	    $attribute = '.';
	    $score = 0;
	}

	print join("\t",
		   $seqname,
		   'avg',
		   $context,
		   $i,
		   $i + $width - 1,
		   sprintf('%.3f', $score),
		   '.',
		   '.',
		   $attribute,
	    ), "\n";
	$i += $step - 1;
    }
    return 0;
}


sub gff_filter_by_coord {
    my ($lower_bound, $upper_bound, $last_index_seen, $data_ref) = @_;

    my @filtered;
    for (my $i = $last_index_seen; $i < @{$data_ref}; $i++) {

        my $start_coord = (split /\t/, $data_ref->[$i])[3];

	if ($start_coord >= $lower_bound && $start_coord <= $upper_bound) {
	    push @filtered, $data_ref->[$i];
	    $last_index_seen = $i;
	}

	last if ( $start_coord > $upper_bound );

    }
    unshift @filtered, $last_index_seen;
    return \@filtered;
}

sub gff_sort {
    my $data_ref = shift;

    print STDERR "sorting data...\n";

    return sort {
	(split /\t/, $a)[0] cmp (split /\t/, $b)[0] or
        (split /\t/, $a)[2] cmp (split /\t/, $b)[2] or
        (split /\t/, $a)[3] <=> (split /\t/, $b)[3]
    } @{$data_ref};
}

sub gff_read {
    my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, shift);
    my %rec = (
	'seqname'=>$seqname,
	'source'=>$source,
	'feature'=>$feature,
	'start'=>$start,
	'end'=>$end,
	'score'=>$score,
	'strand'=>$strand,
	'frame'=>$strand,
	'attribute'=>$attribute
	);
    return \%rec;
}

sub gff_find_array {
    my ($field, $array_ref) = @_;

    print STDERR "finding indices...\n";

    # @index contains a list of locations where array should be split
    # $previousid and $currentid are scalars containing the previous and current IDs (chr1, 2, etc)
    # $chrcount is just a counter for the number of different types of IDs
    my (@index, $previous, $current);
    my $chrcount = 0;

    for (my $i = 0; $i < $#{$array_ref} ; $i++) {

	# gets current starting coordinate
	$current = (split /\t/, $array_ref->[$i])[$field]; #chr1, chr2, etc

	# if we're at beginning of file it doesn't make sense to look for changes already
	# gets previous starting coordinate
	if ($i != 0) {$previous = (split /\t/, $array_ref->[$i-1])[$field];}
	else {$previous = $current;}

	# keeps track of number of different types of records
	# also stores each record type change in @index
	# ignores pound (#) characters if they're the first printing character in the record
	if ( ($current ne $previous) or ($i == 0) )  {
	    $index[$chrcount] = $i;
	    $chrcount++;
	}
    }
    return @index;
}

sub gff_print_header {
    my @call_args = @_;
    print '##gff-version 3\n';
    print join(' ',
	       '#',
	       @call_args,
	       "\n");
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime (time);
    printf "# %4d-%02d-%02d %02d:%02d:%02d\n", $year+1900, $mon+1, $mday, $hour, $min, $sec;
    print join("\t",
	       '# SEQNAME',
	       'SOURCE',
	       'FEATURE',
	       'START',
	       'END',
	       'SCORE',
	       'STRAND',
	       'FRAME',
	       'ATTRIBUTES',
	), "\n";
    return 0;
}


__END__

=head1 NAME

 window_gff.pl - Generate windows from GFF attributes

=head1 VERSION

 $Rev::                                                     $: Revision of last commit
 $Author::                                                  $: Author of last commit
 $Date::                                                    $: Date of last commit
 $HeadURL::                                                 $: URL to last version in repository
 $Id$:

=head1 USAGE

 # window file using implicit standard input and unix pipes and redirection
 cat foo.gff | window_gff.pl --width 200 --step 100 > windowed_bar.gff

 # window file using explicit standard input and short options
 cat bar.gff | window_gff.pl --gff-file - -w 400 -s 200 -o windowed_bar.gff

 # window single file with no overlaps
 window_gff.pl -f foo.gff -w 100 --output windowed_foo.gff

 # window multiple pre-sorted files in batch mode
 window_gff.pl --batch foo.gff bar.gff --width 500 --no-sort

 # window multiple files in batch mode using shell expansion with extra warnings
 window_gff.pl -b *.gff -w 100 --verbose

=head1 REQUIRED ARGUMENTS

 -w, --width    width size of the sliding window in base pairs

=head1 OPTIONS

 countMethylation.pl [OPTION]... [FILE]...

 -f, --gff-file    GFF alignment input file
 -w, --width       width size of sliding window in bp
 -s, --step        step interval of sliding window in bp
 -n, --no-sort     assumes input gff file is pre-sorted by sequence, feature, and coordinate
 -b, --batch       takes any number of filenames as arguments and windows them in batch mode
 -o, --output      filename to write results to (default is STDOUT, unless in batch mode)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 DESCRIPTION

 Takes a GFF-formatted input file with methylation information on the attribute field
 in the form c=? and t=?. Each line corresponds to a single 'c' in a sequenced genome.
 The input file may contain multiple sequence id's and features-contexts.
 Runs a sliding window of width x and step interval y and averages the score
 for each step, generating a GFF file with N/y lines (N=number of input 'c's).
 The averaging is done by extracting the attribute fields and computing score = c / (c + t)

=head1 SUBROUTINES

=head2 gff_sliding_window()

 Takes as input the window width, step interval, and a gff array
 Groups input data into @data/step windows of overlap width/step
 Prints each window as a GFF line

=head2 gff_filter_by_coord()

 Takes min and max coords allowed, the coord to search for, a starting offset
 and an gff array reference.
 Returns reference to array with all lines having coordinates in given range

=head2 gff_sort()

 Takes reference to array to be sorted
 Sorts gff lines by sequence, feature and start coordinate
 Returns sorted array

=head2 gff_read()

 Takes a tab-delimited gff line string
 Returns reference to hash with keys as defined in GFF spec v3

=head2 gff_find_array()

 Takes a field specifier and reference to array to split
 Finds each index that signifies a change in type of record (according to field)
 Returns index array

=head2 gff_print_header()

 Prints a commented line with header fields based on the GFF v3 spec

=head1 DIAGNOSTICS

=head1 CONFIGURATION AND ENVIRONMENT

=head1 DEPENDENCIES

=head1 INCOMPATIBILITIES

=head1 BUGS AND LIMITATIONS

=head1 AUTHOR

 Pedro Silva <psilva@nature.berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 LICENSE AND COPYRIGHT

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
