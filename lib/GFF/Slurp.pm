package GFF::Slurp;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use GFF::Parser;
use GFF::Util;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(gff_slurp gff_slurp_index);


=head2 gff_slurp('file')

return an arrayref of all gff records

=cut

sub gff_slurp{
    my ($file_or_handle) = @_;
    my $p = GFF::Parser->new(file => $file_or_handle);
    my @accum;
    while (my $gff = $p->next()){
        # if we're skipping, no need to check
        if (is_gff($gff)){ 
            push @accum, $gff;
        }
    }
    return \@accum;
}

=head2 gff_slurp_index('file','chr1')

return a hashref of column val to gff record

=cut

sub gff_slurp_index{
    my ($file_or_handle, $column) = @_;
    my $p = GFF::Parser->new(file => $file_or_handle);

    my %index; # { index => [gffs] }

    my $counter = 0;
    my $badrecords = 0;

    while (my $gff = $p->next()){
        my $i = $gff->get_column($column);
        if (defined $i){
            push @{$index{$i}}, $gff;
        } else {
            $badrecords++;
        }
        $counter++;
    }

    for my $i (keys %index) {
        my @gffs = @{$index{$i}};
        $index{$i} = [
        map  { $_->[1] }
        sort { $a->[0] <=> $b->[0] }
        map  { [$_->start,$_] } 
        @gffs
        ];
    }

    carp "warning: $badrecords out of $counter records didn't have a column/attribute $column in " .  $p->filename_or_handle
    if $badrecords;

    return \%index;
}

1;

