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
use Window;
use WindowSource;

Log::Log4perl->easy_init({ 
        level    => $DEBUG,
        layout   => "%d{HH:mm:ss} %p> (%L) %M - %m%n", 
        #file     => ">>log4perl.log",
    });

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ! ($opt_target xor $opt_fixed);

if ($opt_output ne '-'){
    close STDOUT;
    open STDOUT, '>', $opt_output;
}

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
my @sequences = keys %split_queries;

INFO("%split_queries: " . Dumper \%split_queries);

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

        window_gff($query_parser, $target_window_source);
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

        window_gff($query_parser, $target_window_source);
    }
}

sub window_gff{
    my ($query_parser, $target_window_source)  = @_;
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
                flush($win);
            } else {
                push @keep, $i;
            }
        }

        @window_pool = @window_pool[@keep];
    }
    foreach my $win (@window_pool) {
        flush($win);
    }
    if ($opt_no_skip){
        while (my $win = $target_window_source->next()){
            flush($win);
        }
    }
}

sub flush{
    my $win = shift;
    if ($win->counter || $opt_no_skip){
        my $s = $win->to_gff;
        DEBUG("This is being flushed:\n" . $s);
        say $s;
    }
}


=head1 NAME

window_gff_new.pl - Your program here

=head1 SYNOPSIS

Usage examples:

 window_gff_new.pl [options]...

=head1 REQUIRED ARGUMENTS

=for comment #####################################################################

=over

=item --query <file>

The GFF file we are windowing. 

=for Euclid
    file.type:  readable

=item --scoring <method>

Scoring Method to use. (See SCORING METHODS section below for detail).

=back

=head1 OPTIONS

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

=cut

