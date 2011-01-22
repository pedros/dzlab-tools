#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use Pod::Usage;
use List::Util 'shuffle';

my $DATA_HANDLE = 'ARGV';
my $output;
my $loci_list;
my $include       = 0;
my $attribute_id  = 'ID';
my $random_sample = 0;
my $type;
my $min_score     = 0;
my $max_score     = 0;
my $min_length    = 0;
my $max_length    = 0;

# Grabs and parses command line options
my $result = GetOptions(
    'list|l=s'          => \$loci_list,
    'include|i'         => \$include,
    'attribute-id|a=s'  => \$attribute_id,
    'random-sample|r=i' => \$random_sample,
    'min-score|s=f'     => \$min_score,
    'max-score|S=f'     => \$max_score,
    'min-length|m=f'    => \$min_length,
    'max-length|M=f'    => \$max_length,
    'type|t=s'          => \$type,
    'output|o=s'        => \$output,
    'verbose'         => sub { use diagnostics; },
    'quiet'           => sub { no warnings; },
    'help'            => sub { pod2usage( -verbose => 1 ); },
    'manual'          => sub { pod2usage( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage( -verbose => 1 )
    unless @ARGV and $result;

if ($output) {
    open my $USER_OUT, '>', $output
        or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my @random;
my %loci_to_delete;

if ($loci_list) {
    open my $LOCI, '<', $loci_list or croak "Can't read $loci_list: $!";
    while (defined(my $locus = <$LOCI>)){
        if ($locus =~ /(\w+)/){
            $loci_to_delete{$1} = 1;
        }
    }
    close $LOCI or die "Can't close $loci_list: $!";
}

my $gff_iterator = make_gff_iterator( $ARGV[0], \&gff_read );
GFF:
while ( my $gff_line = $gff_iterator->() ) {

    next GFF unless ref $gff_line eq 'HASH';

    my $attribute = $gff_line->{attribute};
    $gff_line->{attribute} =~ s/.*
                                $attribute_id
                                [=\s]?
                                ([^;\s]+)
                                .*
                               /$1/mx
                               if $attribute_id;
    
    my $matches 
    =   (( 0 == $min_score  or $gff_line->{score} eq q{.} or $gff_line->{score} >= $min_score )
    and ( 0 == $max_score  or $gff_line->{score} eq q{.} or $gff_line->{score} <= $max_score )
    and ( 0 == $min_length or ( $gff_line->{end} - $gff_line->{start} ) >= $min_length )
    and ( 0 == $max_length or ( $gff_line->{end} - $gff_line->{start} ) <= $max_length )
    and ( not defined $type or $gff_line->{attribute} =~ m/$type/ )
    and ( not defined $loci_list or exists $loci_to_delete{ $gff_line->{attribute} } ));

    next GFF if $include xor $matches;

    my $out_string = join( "\t",
                           $gff_line->{seqname}, $gff_line->{source}, $gff_line->{feature},
                           $gff_line->{start},   $gff_line->{end},    $gff_line->{score},
                           $gff_line->{strand},  $gff_line->{frame},  $attribute ) . "\n";

    if ($random_sample) {
        push @random, $out_string;
    }
    else {
        print $out_string;
    }
}

if ($random_sample) {
    croak "Not enough loci for sample" if $random_sample > @random;
    @random = shuffle @random;
    print @random[ 0 .. $random_sample - 1 ];
}

sub make_gff_iterator {
    my ( $gff_file_name, $gff_parser ) = @_;

    open my $GFF_FILE, '<', $gff_file_name
        or croak "Can't read $gff_file_name: $!";

    return sub { $gff_parser->( scalar <$GFF_FILE> ) };
}

sub gff_read {
    return [] if $_[0] =~ m/^
                            \s*
                            \#+
                           /mx;

    my ($seqname, $source, $feature, $start, $end,
        $score,   $strand, $frame,   $attribute
    ) = split m/\t/xm, shift || return;

    $attribute =~ s/[\r\n]//mxg;

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

__END__


=head1 NAME

 filter_gff.pl - Short description

=head1 SYNOPSIS

 filter_gff.pl -l loci_to_delete.txt -a ID annotation_to_filter.gff

=head1 DESCRIPTION

=head1 OPTIONS

 filter_gff.pl [OPTION]... [FILE]...

 -l, --list          list of loci IDs, one ID per line
 -s, --min-score     minimum score to filter by (0 by default)
 -S, --max-score     maximum score to filter by (0 by default)
 -m, --min-length    minimum length to filter by (0 by default)
 -M, --max-length    maximum length to filter by (0 by default)
 -i, --include       only print out those loci in GFF with matching IDs to the loci-list
 -a, --attribute-id  attribute tag in GFF file (default 'ID')
 -t, --type          match type in attribute value
 -r, --random-sample take a random sample of size of the parameter
 -o, --output        filename to write results to (defaults to STDOUT)
 -v, --verbose       output perl's diagnostic and warning messages
 -q, --quiet         supress perl's diagnostic and warning messages
 -h, --help          print this information
 -m, --manual        print the plain old documentation page

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
