#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;

pod2usage(-verbose => 99,-sections => [qw/NAME SYNOPSIS OPTIONS/]) 
if ! ($opt_target xor $opt_fixed);

use FindBin;
use lib "$FindBin::Bin/lib";
use GFF::Parser::Attributes;

my $query_parser = GFF::Parser::Attributes->new(file => $opt_query);

=head1 NAME

window_gff_new.pl - Your program here

=head1 SYNOPSIS

Usage examples:

 window_gff_new.pl [options]...

=head1 REQUIRED ARGUMENTS

=over

=item --query <file>

The thing we are windowing. Single sequnce gff sorted by start column.

=for Euclid
    file.type:  readable

=back

=head1 OPTIONS

Fixed Window options:

=over

=item  --fixed <window_size>

Size of fixed window. Cannot be combined with --target

=for Euclid
    window_size.type:    int >= 1 

=item  --step <step_size>

Step size between window (ie, window1.start - window2.start). Defaults to window size 


=for Euclid
    step_size.type:    int >= 1 

=back

Annotation options:

=over

=item --target <file>

Single sequnce gff sorted by start column. Cannot be combined with --fixed

=for Euclid
    file.type:  readable

=item --output <file>

Output

=for Euclid
    file.type:  string
    file.default:  '-'

=item --genome <fasta>

Original reference genome.  This is used only for acquiring lengths of the chromosomes.

=for Euclid
    file.type:  readable

=item --help

=back

=cut

