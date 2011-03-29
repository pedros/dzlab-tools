#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use autodie;
use Getopt::Euclid qw( :vars<opt_> );
use Pod::Usage;
use feature 'say';

pod2usage(-verbose => 99,-sections => [('NAME', 'SYNOPSIS', 'OPTIONS', 'REQUIRED ARGUMENTS')] )
    if $opt_help;

open my $ain, '<', $opt_input_a;
open my $bin, '<', $opt_input_b;

open my $aout, '>', $opt_output_a;
open my $bout, '>', $opt_output_b;

open my $aerror, '>', "$opt_output_a.$opt_error_suffix";
open my $berror, '>', "$opt_output_b.$opt_error_suffix";

CMP:
while (defined (my $a_record = <$ain>) and 
    defined (my $b_record = <$bin>)) {

    my ($a_id, $a_mm, $a_rawcoord) = (split /\t/, $a_record)[0,2,3];
    my ($b_id, $b_mm, $b_rawcoord) = (split /\t/, $b_record)[0,2,3];

    $a_rawcoord =~ s/.*chr\w:(\d+).*/$1/xmsi;
    $b_rawcoord =~ s/.*chr\w:(\d+).*/$1/ixms;

    $a_mm = get_score ($a_mm);
    $b_mm = get_score ($b_mm);

    #say $a_rawcoord;
    #say $b_rawcoord;

    if (! defined $a_mm and ! defined $b_mm) {
        next CMP; # no matches at all
    }
    elsif (defined $a_mm and defined $b_mm){
        next CMP if $a_mm == $b_mm;
        if ($a_rawcoord == $b_rawcoord){
            print $aout $a_record if $a_mm < $b_mm;
            print $bout $b_record if $a_mm > $b_mm;
        } else {
            print $aerror $a_record;
            print $berror $b_record;
        }
    }
    elsif (! defined $a_mm) {
        print $bout $b_record;
    }	
    elsif (! defined $b_mm) {
        print $aout $a_record;
    }
    else {croak "Impossible situation:\n$a_record\n$b_record"}
}

close $aout; close $bout;
close $ain;   close $bin;
close $aerror;
close $berror;

sub get_score {
    my ($mm) = @_;

    return if 'NM' eq $mm;
    my @mm = split /:/, $mm;

    for my $i (0 .. @mm - 1) {
        return $i if 1 == $mm[$i];
    }
}

__END__


=head1 NAME

 split_on_mismatches_2.pl - Filter sorted bowtie inputs into two files based on mismatch counts

=head1 SYNOPSIS

 # overwrites input files with correct imprinting alignments
 perl split_on_mismatches_2.pl 
    -ea Col -eb Ler 
    -a CxL_En_WT_Col.bowtie -b CxL_En_WT_Ler.bowtie 
    -oa filtered_a -ob filtered_b

=head1 REQUIRED ARGUMENTS

=over

=item  -ia <eland> | --input-a <eland>

Ecotype A eland file

=for Euclid eland.type:        readable

=item  -ib <eland> | --input-b <eland>

Ecotype B eland file

=for Euclid eland.type:        readable

=item  -oa <output> | --output-a <output>

Filtered A eland output file

=item  -ob <output> | --output-b <output>

Filtered B eland output file

=back

=head1 OPTIONS

=over

=item --coord | -c

Use coordinates for matching as well.

=item  -e <suffix> | --error-suffix <suffix>

=for Euclid
    suffix.default:     'error'


=item --help | -h

=back

=cut

