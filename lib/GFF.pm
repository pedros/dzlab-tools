package GFF;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Moose;
use Carp;
use autodie;    

has asterisk         => ( is => 'rw', isa => 'Bool', default => 0);
# type constraints very slow? removed.
has sequence         => ( is => 'ro');
has source           => ( is => 'ro');
has feature          => ( is => 'ro');
has start            => ( is => 'ro');
has end              => ( is => 'ro');
has score            => ( is => 'ro');
has strand           => ( is => 'ro');
has frame            => ( is => 'ro');
has attribute_string => ( is => 'rw');
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
    my $attrstr = $self->attribute_string;
    if (defined($attrstr) ){
        if ($attrstr=~s/\*$//){
            $self->asterisk(1);
            $self->attribute_string($attrstr);
        }
        for (split /;/, $attrstr){
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
