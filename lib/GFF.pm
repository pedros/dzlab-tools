package GFF;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    

has sequence    => ( is => 'ro', isa => 'Maybe[Str]',);
has source      => ( is => 'ro', isa => 'Maybe[Str]',);
has feature     => ( is => 'ro', isa => 'Maybe[Str]',);
has start       => ( is => 'ro', isa => 'Maybe[Int]',);
has end         => ( is => 'ro', isa => 'Maybe[Int]',);
has score       => ( is => 'ro', isa => 'Maybe[Str]',);
has strand      => ( is => 'ro', isa => 'Maybe[Str]',);
has frame       => ( is => 'ro', isa => 'Maybe[Str]',);
has attribute_string => ( is => 'ro', isa => 'Maybe[Str]',);
has attr_hash => (
      traits    => ['Hash'],
      is        => 'ro',
      isa       => 'HashRef[Str]',
      handles   => {
          set_attribute  => 'set',
          get_attribute  => 'get',
          list_attribute => 'keys',
      },
      lazy_build => 1,
  );

sub _build_attr_hash{
    my $self = shift;
    my %accum;
    if (defined($self->attribute_string) ){
        for (split /;/, $self->attribute_string){
            my ($key, $val) = split /=/, $_;
            $key =~ s/^\s+//;
            $key =~ s/\s+$//;

            if (defined $val){
                $val =~ s/^\s+//;
                $val =~ s/\s+$//;
                $accum{$key} = $val;
            }
            else {
                $accum{Note} = $key;
            }
        }
    }
    return \%accum;
}

sub to_string{
    my $self = shift;
    return 
    join "\t",
    map { ! defined $_ ? q{.} : $_ } 
    ($self->sequence, $self->source, $self->feature, $self->start, $self->end,
        $self->score,   $self->strand, $self->frame,   $self->attribute_string);
}

sub parse_locus{
    my ($self, $locus_tag) = @_;
    if ($self->get_attribute($locus_tag) =~ /([^\.]+)\.([^\.]+)/){
        return ($1,$2);
    } else{
        return;
    }
}

my %cols = map { $_ => 1 } qw/sequence source feature start end score strand frame attribute_string/;

sub get_column{
    my ($self, $colname) = @_;
    if (exists $cols{$colname}){
        return $self->$colname;
    } else{
        return $self->get_attribute($colname);
    }
}


sub equals{
    my ($self, %against) = @_;
    my @columns = keys %against;
    foreach my $col (@columns) {
        my $x = $self->get_column($col);
        my $y = $against{$col};

        if ((!defined $x && !defined $y) || (defined $x && defined $y && $x eq $y)){
            next;
        } else {
            return 0;
        }
    }
    return 1;
}

sub length{
    my $self = shift;
    return $self->end - $self->start + 1;
}

no Moose;
__PACKAGE__->meta->make_immutable;

#sub colname2num{ return $colmap{$_[0]} // croak "bad colmn name"; }

1;
