package DZLab::Tools::GFF;

=head1 NAME

DZLab::Tools::GFF - IO and parsing utilities for GFFv3

=head1 VERSION

This document describes DZLab::Tools::GFF version 0.0.1

=head1 SYNOPSIS

    use DZLab::Tools::GFF qw/gff_read gff_make_iterator/;

    my $iterator = gff_make_iterator( parser => \&gff_read, handle => \*STDIN );

    # or your own parser:
    $iterator = gff_make_iterator( sub { split /\t/ }, file => 'file.gff');

    while (my $gff = iterator->()) {
        # do stuff with $gff
    }

=head1 DESCRIPTION

Using BioPerl for simple parsing of GFF3 files can be overkill, given its dependencies. 
This module offers a functional approach to the problem.

=cut

use strict; use warnings;
use version; our $VERSION = '0.0.1';
use Carp;
require Exporter;

our @ISA       = qw/Exporter/;
our @EXPORT    = qw//;
our @EXPORT_OK = qw/gff_read gff_make_iterator/;

=head1 EXPORTED FUNCTIONS

=head2 gff_read $string

Parses a single GFF3 line.

Returns undef or empty list on undefined $string (EOF).

Returns [] on blank lines and comments.

Returns [pragma, value, ...] on comments of form '##pragma value value ...'

Returns reference to hash with keys otherwise:

    seqname       => string | .
    source        => string | .
    feature       => string | .
    start         => int
    end           => int
    score         => float  | .
    strand        => + | - | .
    frame         => 0 | 1 | 2 | .
    attribute     => string
    attributes    => {key => string | [string, ...], ...}

=cut
sub gff_read {
    my ($gff_line) = @_;
    return unless defined $gff_line;

    # get pragmas of the form: '##gff-version 3' and '##sequence-region ctg123 1 1497228'
    return [split /\s+/, $1] if $gff_line =~ m/^
                                               \s*
                                               \#{2}
                                               \s*
                                               (.*)
                                               $/mx;

    # ignore blank lines and lines starting with '#'
    return [] if $gff_line =~ m/^ \s* (?:\#+ .*)? $/mx;

    my ($seqname, $source, $feature, $start, $end,
        $score,   $strand, $frame,   $attribute
    ) = split m/\t/xm, $gff_line || return;

    $attribute =~ s/[\r\n]//mxg;

    my %attributes = map { /=/ ? split /=/ : (Note => $_) } split /;/, $attribute;
    @attributes{keys %attributes} = map { /,/ ? [split /,/] : $_ } values %attributes;

    return {
        seqname   => lc $seqname,
        source    => $source,
        feature   => $feature,
        start     => $start,
        end       => $end,
        score     => $score,
        strand    => $strand,
        frame     => $frame,
        attribute => $attribute,
        attributes=> \%attributes
    };
}

=head2 gff_make_iterator %options

Returns an anonymous function that on, each call, reads one gff line from a $file or $handle, and returns a parsed structure as delivered by $parser.

Returns undef or empty list on EOF.

=cut
sub gff_make_iterator {
    my %options = @_;

    my $parser = $options{parser};
    my $file   = $options{file};
    my $handle = $options{handle};

    croak
    "Need parser function reference and file name or handle to build iterator"
    unless $parser and ref $parser eq 'CODE'
    and (
        (defined $file and -e $file)
        or
        (defined $handle and ref $handle eq 'GLOB')
    );

    if (defined $file and -e $file) {
        undef $handle;
        open $handle, '<', $file
        or croak "Can't read $file: $!";
    }

    return sub {
        $parser->( scalar <$handle> );
    };
}


1;

=head1 INSTALLATION

To install this module type the following:

   perl Build.PL
   Build
   Build test
   Build install

or

   perl Makefile.PL
   make
   make test
   make install

=head1 DIAGNOSTICS

=over

=item C<< "Need parser function reference and file name or handle to build iterator" >>

C<< gff_make_iterator >> takes a hash of options. 'parser', and either 'file' of 'handle' is mandatory.
'parser' must be a CODE ref, file must be defined and exist in the filesystem, and handle must be a GLOB ref.

=back

=head1 TODO

=over

=item

Add a GFF3 validator

=item

Add a GFF3 pretty-printer

=item

Add an output formatting function

=item

Merge other GFF-using subroutines into this module

=back

=head1 AUTHOR

Pedro Silva  C<< <pedros@berkeley.edu> >>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010, Pedro Silva C<< <pedros@berkeley.edu> >>. All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see L<http://www.gnu.org/licenses/>.


