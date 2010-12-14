package DZLab::Tools::BinaryRangeSearch;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(binary_range_search);

use DZLab::Tools::RangeUtils;

=head1 binary_range_search {queries => [...ranges...], database =>
[...gff_records...]}

given a list of query ranges, and a list of database gff_records, each with
it's own range, return an iterator which iterates over gff_records ranges 
that overlaps with a query range

=cut

sub binary_range_search {
    my %options = @_;

    my $queries  = $options{queries}  || croak 'Need a queries parameter';
    my $database = $options{database} || croak 'Need a database parameter';

    my ( $low, $high ) = ( 0, $#{$database} );
    my @iterators = ();

TARGET:
    for my $query (@$queries) {

    RANGE_CHECK:
        while ( $low <= $high ) {

            # middle
            my $try = int( ( $low + $high ) / 2 );

            if (range_before($database->[$try]{range}, $query)){
                $low = $try + 1;
                next RANGE_CHECK;
            } 
            elsif (range_before($query, $database->[$try]{range})){
                $high = $try - 1;
                next RANGE_CHECK;
            }

            # if we make it here, we've found a source $try which overlaps with the  
            # target $query
            
            my ( $down, $up ) = ($try) x 2;
            my %seen = ();

            # create an iterator which, on every call, returns an overlapping
            # source range, starting with $query and then going down or up the
            # list.
            my $brs_iterator = sub {
                
                # if $query overlaps with $up + 1, and it's new, return it
                if ( range_overlap($database->[ $up + 1 ]{range}, $query)
                    # $database->[ $up + 1 ]{end} >= $query->[0]
                    # and $database->[ $up + 1 ]{start} <= $query->[1]
                    and !exists $seen{ $up + 1 } )
                {
                    $seen{ $up + 1 } = undef;
                    return $database->[ ++$up ];
                }
                # if $query overlaps with $down - 1, and it's new, return it
                elsif ( range_overlap($database->[$down-1]{range},$query)
                    #$database->[ $down - 1 ]{end} >= $query->[0]
                    #and $database->[ $down - 1 ]{start} <= $query->[1]
                    and !exists $seen{ $down - 1 }
                    and $down > 0 )
                {
                    $seen{ $down - 1 } = undef;
                    return $database->[ --$down ];
                }
                # we already know $try overlaps, so return it too.
                elsif ( !exists $seen{$try} ) {
                    $seen{$try} = undef;
                    return $database->[$try];
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
        while (@iterators) {
            if ( my $range = $iterators[0]->() ) {
                return $range;
            }
            shift @iterators;
        }
        return;
        };
}

1;

