package GFF;
use Data::Dumper;
use feature 'say';

use warnings;
use strict;
use Carp;

use FindBin;
use lib "$FindBin::Bin";
use GFF::Parser;

use version; our $VERSION = qv('0.0.1');

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(colname2num);
our @EXPORT = qw(gff_to_string);

=head1 EXPORT

=head1 SUBROUTINES/METHODS

=head2 gff_to_string $gffrec

return a gffrec back in original text format

=cut

sub gff_to_string{
    my ($gff) = @_;
    croak "gff_to_string needs an argument??" unless $gff;
    if (ref $gff eq 'HASH'){
        return 
        join "\t",
        map { ! defined $_ ? q{.} : $_ } 
            @{$gff}{'seqname', 'source', 'feature', 'start', 'end',
            'score',   'strand', 'frame',   'attribute'};
    } 
    else {
        croak "non-gff record given to gff_to_string";
    }
}
my %colmap = (
    seqname   => 1,
    source    => 2,
    feature   => 3,
    start     => 4,
    end       => 5,
    score     => 6,
    strand    => 7,
    frame     => 8,
    attribute => 9
);

sub colname2num{ return $colmap{$_[0]} // croak "bad colmn name"; }

=head1 AUTHOR

tnish, C<< <tnish at berkeley.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-gff at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=GFF>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc GFF

=head1 ACKNOWLEDGEMENTS

=head1 LICENSE AND COPYRIGHT

Copyright 2011 tnish.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of GFF
