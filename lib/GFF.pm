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
our @EXPORT_OK = qw();
our @EXPORT = qw(is_gff gff_to_string);

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
    elsif (ref $gff eq 'ARRAY'){
        return 
        join "\t",
        map { ! defined $_ ? q{.} : $_ } 
        @{$gff}[0..8];
    }
    else {
        croak "non-gff record given to gff_to_string";
    }
}

=head2 do_gff BLOCK ARGS_FOR_GFFPARSER_CONSTRUCTOR

BLOCK is executed with $_ set to the GFF hashref.  ARGS_FOR_GFFPARSER_CONSTRUCTOR are passed to
GFF::Parser constructor directly.

 do_gff { 
    # ... $_ is a gff hashref here
 } file => 'file.gff', locus => 'ID';

=cut

sub do_gff(&@){
    my ($code, @opt) = @_;
    my $iter = GFF::Parser->new(@opt);
    while (my $gff = $iter->next()){
        local $_ = $gff;
        &$code;
    }
}

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
