#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use Log::Log4perl qw/:easy/;
use FindBin;
use lib "$FindBin::Bin/lib";
use Fasta;
use GFF::Parser;
use GFF::Split;
use GFF::Sort;
use Window;
use WindowSource;

Log::Log4perl->easy_init({ 
        level    => $DEBUG,
        layout   => "%d{HH:mm:ss} %p> (%L) %M - %m%n", 
        file     => ">log4perl.log",
    });

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS SCORING/]) 
if ! ($opt_target && $opt_locus_tag xor $opt_fixed && $opt_genome);


### Scoring

my %scoring_dispatch = (
    sum             => 'Window::SumScore',
    fracmeth        => 'Window::FractionalMethylation',
    average         => 'Window::AverageScore::Full',
    diluted_average => 'Window::AverageScore::Diluted',
    locus           => 'Window::Locus',
);
my $scoring = $scoring_dispatch{$opt_scoring};
INFO("Scoring Method: $scoring");

### Split 

INFO("Splitting query file by sequence and sorting each by start position.");

my %split_queries = gff_split(file => $opt_query, sort => 'start', sequence => 'all', keep => 0);
my %temp_outfiles   = map {$_ => $split_queries{$_} . '.windowed'} keys %split_queries;
my %temp_outfh      = map {open my $fh, '>', $temp_outfiles{$_}; $_ => $fh} keys %temp_outfiles;

my @sequences = keys %split_queries;

INFO("%split_queries: " . Dumper \%split_queries);
INFO("%temp_output: "   . Dumper \%temp_outfiles);
INFO("%temp_fh: "       . Dumper \%temp_outfh);

INFO("Ok, windowing " . join ", ", @sequences);

if ($opt_target){
    INFO("Starting annotation windowing");

    my %split_target = gff_split(file => $opt_target, sort => 'start', sequence => 'all', keep => 0);

    INFO("%split_target: " . Dumper \%split_target);

    for my $sequence (@sequences){
        my $query_file = $split_queries{$sequence};
        my $target_file = $split_target{$sequence};
        INFO("Windowing $sequence ($query_file, $target_file)");

        INFO("Creating query_parser and window source");
        my $query_parser = GFF::Parser->new(file => $query_file);
        my $target_window_source = WindowSource::Annotation->new(
            class => $scoring,
            file => $target_file,
            locus => $opt_locus_tag,
            sequence => $sequence,
        );

        window_gff($query_parser, $target_window_source, $temp_outfh{$sequence});
    }
}
if ($opt_fixed){
    INFO("Starting fixed windowing");
    INFO("Getting lengths of chromosomes from reference genome");

    my %seq_lengths = %{count_fasta($opt_genome)};
    INFO("%seq_lengths: " . Dumper \%seq_lengths);

    for my $sequence (@sequences){
        my $query_file = $split_queries{$sequence};
        INFO("Windowing $sequence ($query_file)");

        INFO("Creating query_parser and window source");
        my $query_parser = GFF::Parser->new(file => $query_file);
        my $target_window_source = WindowSource::Fixed->new(
            class => $scoring,
            sequence => $sequence,
            length => $seq_lengths{$sequence},
            step => $opt_fixed,
        );

        window_gff($query_parser, $target_window_source, $temp_outfh{$sequence});
    }
}

while (my ($seq,$fh) = each %temp_outfh) {
    INFO("closing $seq file handle");
    close $fh;
}

my $outfh;
if ($opt_output eq '-'){
    $outfh = \*STDOUT;
} else {
    open $outfh, '>', $opt_output;
}

for my $seq (@sequences){
    my $file = $temp_outfiles{$seq};
    INFO("sorting and writing $seq ($file) to $opt_output");
    gff_sort(file => $file, overwrite => 1, column => 'start');

    open my $fh, '<', $file;
    while (defined (my $line = <$fh>)) { print $outfh $line; }
    close $fh;
}
if ($opt_output ne '-'){
    close $outfh;
}

### Done

sub window_gff{
    my ($query_parser, $target_window_source,$fh)  = @_;
    my @window_pool;

    while (my $query = $query_parser->next()){
        DEBUG("Pulled from query_parser: \n" . $query->to_string);
        while (my $w = $target_window_source->next()){
            push @window_pool, $w;
            last if $w->start > $query->end; 
            # (we've already gotten the last possible window $query could overlap with)
        }
        DEBUG("Pool has " . scalar @window_pool . " elements:\n" . join("\n", map {$_->to_gff} @window_pool));

        my @keep; # indices to keep
        for my $i (0 .. $#window_pool){
            my $win = $window_pool[$i];
            if ($win->overlap($query)){
                $win->accumulate($query);
                push @keep, $i;
            } elsif ($win->end < $query->start){
                # pool's end is behind the queries start, so no more will ever match. flush
                flush($win,$fh);
            } else {
                push @keep, $i;
            }
        }

        @window_pool = @window_pool[@keep];
    }
    foreach my $win (@window_pool) {
        flush($win,$fh);
    }
    if ($opt_no_skip){
        while (my $win = $target_window_source->next()){
            flush($win,$fh);
        }
    }
}

sub flush{
    my ($win,$fh) = @_;
    if ($win->counter || $opt_no_skip){
        my $s = $win->to_gff;
        DEBUG("This is being flushed:\n" . $s);
        say $fh $s;
    }
}


=head1 NAME

window_gff_new.pl - Your program here

=head1 SYNOPSIS

 perl window_gff_new.pl --query tiny-query.gff --scoring sum --fixed 50 --genome tiny-genome.fas --output output.txt

=head1 OPTIONS

=for comment #####################################################################

Required options:

=over

=item --query <file>

The GFF file we are windowing. 

=for Euclid
    file.type:  readable

=item --scoring <method>

Scoring Method to use. (See SCORING METHODS section below for detail).

=back

=for comment #####################################################################

Fixed Window options:

=over

=item  --fixed <window_size>

Size of fixed window. Cannot be combined with --target

=for Euclid
    window_size.type:    int >= 1 

=item --genome <fasta>

Original reference genome.  This is used only for acquiring lengths of the chromosomes.

=for Euclid
    fasta.type:  readable

=back

=for comment #####################################################################

Annotation options:

=over

=item --target <file>

Gff annotation file. Cannot be combined with --fixed

=for Euclid
    file.type:  readable

=item --locus-tag <locus_tag>

Locus tag in annotation file.  Used to label each line of output file.

=back

=for comment #####################################################################

Common Options

=over


=item --output <file>

Output

=for Euclid
    file.type:  string
    file.default:  '-'

=item --no-skip

Pass this option if you want windows without any matches to be printed as well.

=item --help

=back

=head1 SCORING

possible values for --scoring are:

=over

=item sum

For each window, sum the score from each overlap and report in the score column (column 6). 

=item fracmeth

For each window, sum the 'n', 'c', and 't' from each overlap and report the fractional methylation c/(c+t) as the score. 

=item average         

For each window, calculate the mean, standard deviation, variance for the scores of overlaps.

=item diluted_average 

For each window, calculate the mean, standard deviation, variance for the DILUTED scores.  Diluted means that 

 Query:    |-------------------|                 Length: x, Score: n
 Window:                |--------------------|   Length: y
 Overlap:               |------|                 Length: z

Then the score contribution of the query to the window is n * (x/z) * (y/z).  This was yvonne's idea so if it doesn't
make sense, blame her.

=item locus          



=back

=cut

