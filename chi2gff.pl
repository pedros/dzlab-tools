#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/reduce/;

my $output;
my $gff;

# Grabs and parses command line options
my $result = GetOptions (    
    'gff|g=s'     => \$gff,
    'output|o=s'  => \$output,
    'verbose|v'   => sub { use diagnostics; },
    'quiet|q'     => sub { no warnings; },
    'help|h'      => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'    => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my %gff_records = ();
if ($gff) {
    
    my $gff_iterator
    = make_gff_iterator( parser => \&gff_read, file => $gff );

  LOAD:
    while ( my $gff_line = $gff_iterator->() ) {
        next LOAD unless ref $gff_line eq 'HASH';
        push @{ $gff_records{ $gff_line->{seqname} } }, $gff_line;
    }
}

print "##gff-version 3\n";
while (<>) {
    s/[\r\n]//g;;
    my @fields    = split /\t/;
    my $attribute = "ID=$fields[0]; p=$fields[5]; maxwin=$fields[6]";

    if ($gff) {
        my $brs_iterator = binary_range_search(
            range    => [ [$fields[2], $fields[3]] ],
            ranges   => $gff_records{lc $fields[1]},
        );

        my @ranges;
        while (my $range = $brs_iterator->()) {
            push @ranges, $range;
        }

        my $max_window = reduce { $a->{score} > $b->{score} ? $a : $b } @ranges;
        $attribute .= "; maxstart=$max_window->{start}; maxend=$max_window->{end}";
    }

    print join ("\t",
                $fields[1],
                'chipotle',
                'locus',
                $fields[2],
                $fields[3],
                $fields[4],
                q{.},
                q{.},
                $attribute,
            ), "\n";
}


sub binary_range_search {
    my %options = @_;

    my $targets = $options{range}  || croak 'Need a range parameter';
    my $ranges  = $options{ranges} || croak 'Need a ranges parameter';

    my ( $low, $high ) = ( 0, $#{$ranges} );
    my @iterators      = ();

  TARGET:
    for my $range ( @$targets ) {

      RANGE_CHECK:
        while ( $low <= $high ) {

            my $try = int( ( $low + $high ) / 2 );

            $low  = $try + 1, next RANGE_CHECK if $ranges->[$try]{end}   < $range->[0];
            $high = $try - 1, next RANGE_CHECK if $ranges->[$try]{start} > $range->[1];

            my ( $down, $up ) = ($try) x 2;
            my %seen      = ();
        
            my $brs_iterator = sub {

                if (    $ranges->[ $up + 1 ]{end}       >= $range->[0]
                        and $ranges->[ $up + 1 ]{start} <= $range->[1]
                        and !exists $seen{ $up + 1 } ) {
                    $seen{ $up + 1 } = undef;
                    return $ranges->[ ++$up ];
                } 
                elsif ( $ranges->[ $down - 1 ]{end}       >= $range->[0]
                          and $ranges->[ $down + 1 ]{start} <= $range->[1]
                          and !exists $seen{ $down - 1 }
                          and $down > 0 ) {
                    $seen{ $down - 1 } = undef;
                    return $ranges->[ --$down ];
                } 
                elsif ( !exists $seen{$try} ) {
                    $seen{$try} = undef;
                    return $ranges->[$try];
                }
                else {
                    return;
                }
            };
            push @iterators, $brs_iterator;
            next TARGET;
        }
    }

    # In scalar context return master iterator that iterates over the list of range iterators.
    # In list context returns a list of range iterators.
    return wantarray 
    ? @iterators 
    : sub { 
        while( @iterators ) {
            if( my $range = $iterators[0]->() ) {
                return $range;
            }
            shift @iterators;
        }
        return;
    }; 
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

sub make_gff_iterator {
    my %options = @_;

    my $parser     = $options{parser};
    my $file       = $options{file};
    my $GFF_HANDLE = $options{handle};

    croak
        "Need parser function reference and file name or handle to build iterator"
        unless $parser and ref $parser eq 'CODE'
        and ( 
                (defined $file and -e $file)
                xor 
                (defined $GFF_HANDLE 
                     and (ref $GFF_HANDLE eq 'GLOB' or $GFF_HANDLE eq 'ARGV'))
            );

    if ($file) {
        open $GFF_HANDLE, '<', $file
            or croak "Can't read $file: $!";
    }

    return sub {
        $parser->( scalar <$GFF_HANDLE> );
    };
}


__END__


=head1 NAME

 chi2gff.pl -- Converts ChIPotle output files to gff format

=head1 SYNOPSIS

 # find maximum window coordinates in the original data, for each peak
 ./chi2gff.pl chipotle-out.peaks -o chipotle-out.gff -g original_data.gff

 -g is optional

=head1 DESCRIPTION


=head1 OPTIONS

=head1 REVISION

 Version 0.0.2

 $Rev: 261 $:
 $Author: psilva $:
 $Date: 2010-01-14 19:08:08 -0800 (Thu, 14 Jan 2010) $:
 $HeadURL: http://dev.dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/chi2gff.pl $:
 $Id: chi2gff.pl 261 2010-01-15 03:08:08Z psilva $:

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
