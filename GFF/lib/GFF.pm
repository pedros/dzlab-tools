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
our @EXPORT = qw(do_gff is_gff gff_to_string);

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

=head2 gff_slurp

Returns arrayref of gff-records 

=cut

sub gff_slurp{
    my ($file_or_handle, $locus) = @_;
    
    my @accum;
    do_gff {
        push @accum,$_;
    } $file_or_handle, $locus;

    return \@accum;
}

=head2 gff_slurp_by_seq 

Returns hash of sequences to gff-records { sequence => [$gff_hashes] }. 

=cut

sub gff_slurp_by_seq {
    my ($file_or_handle, $locus, $sort) = @_;

    my %gff_records = ();
    do_gff {
        my $seq = $_->{seqname} ;
        push @{ $gff_records{ $seq } }, $_;
    } $file_or_handle, $locus;
    
    if ($sort){
        foreach my $seq (keys %gff_records) {
            @{ $gff_records{$seq} }
            = sort { $a->{start} <=> $b->{start} } @{ $gff_records{$seq} };
        }
    }
    return \%gff_records;
}

sub is_gff{
    my ($gff) = @_;
    return ref $gff && ref $gff eq 'HASH';
}

=head2 do_gff(&@)

 do_gff { 
    # ... $_ is a gff hashref here
 } 'file_or_handle', 'locus_tag';

=cut

sub do_gff(&@){
    my ($code, $file_or_handle, $locus) = @_;
    my $iter = GFF::Parser->new(file => $file_or_handle, locus => $locus);
    while (defined (my $gff = $iter->next())){
        if (is_gff($gff)){
            local $_ = $gff;
            &$code;
        }
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
