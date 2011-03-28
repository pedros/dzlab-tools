package GFF::Util;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use GFF;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(parse_gff is_gff);

=head2 parse_gff($line)

read a line, and return a gff hashref.  maps '.' dot columns and non-existenct attributes to undef. returns 0 on
non-parseable lines or comments. *Returns pragmas as strings*.  Check that ref $result eq 'HASH'

it's a member method so it can be inherited

=cut

sub parse_gff{
    my ($line) = @_;
     
    return 0 unless $line;
    $line =~ tr/\n\r//d;

    if ($line =~ m/^\s*#(#)?/o){
        return defined $1 ? $' : 0; # $' = regex postmatch
    }

    my @split = split /\t/, $line;
    (carp "unparseable GFF line: $line" && return 0) unless @split == 9;

    my %accum;
    @accum{qw/sequence source feature start end score strand frame attribute_string/}
    = map { $_ eq q{.} ? undef : $_ } @split;

    if (defined $accum{sequence}){
        $accum{sequence} = lc $accum{sequence};
    }

    return GFF->new(%accum);
}

sub is_gff{
    my ($gff) = @_;
    return ref $gff eq 'GFF';
}

1;
