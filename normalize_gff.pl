#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use version; our $VERSION = qv('0.0.1');

use FindBin;
use lib "$FindBin::Bin/DZLab-Tools/lib";
use DZLab::Tools::RunUtils;

GetOptions(
    \%ARGV,
    'input|i=s', 'output|o=s', 'error|e=s',
    'gff-a|a=s', 'gff-b|b=s', 'scoring|c=s',
    'feature|f=s', 'dot-as-zero|0',
    _meta_options( \%ARGV ),
) and (@ARGV or $ARGV{input} or @ARGV{qw/gff-a gff-b/}) or pod2usage( -verbose => 1 );

my ( $INH, $OUTH, $ERRH ) = _prepare_io( \%ARGV, \@ARGV );

my %scoring_dispatch = (
    log2_ratio => \&log2_ratio,
    difference => sub { 
        (q{.} eq $_[0] or q{.} eq $_[1])
        ? q{.}
        : ($_[0] - $_[1])
    },
    multiplication => sub {
      (q{.} eq $_[0] or q{.} eq $_[1])
        ? q{.}
        : ($_[0] * $_[1])
    },
);

my $gff_a_iterator = make_gff_iterator( parser => \&gff_read, file => $ARGV{'gff-a'} );
my $gff_b_iterator = make_gff_iterator( parser => \&gff_read, file => $ARGV{'gff-b'} );

GFF:
while (     defined( my $gff_a = $gff_a_iterator->() )
        and defined( my $gff_b = $gff_b_iterator->() ) ) {

    next GFF unless 'HASH' eq ref $gff_a and 'HASH' eq ref $gff_b;

    die "GFF record mismatch: ", Dumper [$gff_a, $gff_b]
    unless @{$gff_a}{qw/seqname start end/} eq @{$gff_b}{qw/seqname start end/};

    if ($ARGV{'dot-as-zero'}) {
        $gff_a->{score} = 0 if $gff_a->{score} eq q{.};
        $gff_b->{score} = 0 if $gff_b->{score} eq q{.};
    }

    print $OUTH join( "\t",
                      $gff_a->{seqname},
                      $gff_a->{source},
                      ($ARGV{feature} // sprintf( "%s/%s", $gff_a->{feature}, $gff_b->{feature} )),
                      $gff_a->{start},
                      $gff_a->{end},
                      (sprintf( "%g" , $scoring_dispatch{$ARGV{scoring}}->( $gff_a->{score}, $gff_b->{score} ))),
                      $gff_a->{strand},
                      $gff_a->{frame},
                      (sprintf( "%s; %s", $gff_a->{attribute}, $gff_b->{attribute})),
                  ), "\n";
}

sub log2_ratio {
    my ($score_a, $score_b) = @_;

    return q{.} if q{.} eq $score_a or q{.} eq $score_b;
    return 'NaN' unless 0 < $score_a and 0 < $score_b;
    return (log( $score_a/$score_b ) / log( 2 ));
}


sub make_gff_iterator {
    my %options = @_;

    my $parser = $options{parser};
    my $file   = $options{file};

    croak
    "Need parser function reference and file name or handle to build iterator"
    unless $parser
    and ref $parser eq 'CODE'
    and (defined $file and -e $file);
    
    open my $GFF_HANDLE, '<', $file
    or croak "Can't read $file: $!";
    
    return sub {
        $parser->( scalar <$GFF_HANDLE> );
    };
}

sub gff_read {
    return [] if $_[0] =~ m/^
                            \s*
                            \#+
                           /mx;

    my ($seqname, $source, $feature, $start, $end,
        $score,   $strand, $frame,   $attribute
    ) = split m/\t/xm, shift || return;

    $attribute =~ s/[\r\n]//mxg;

    return {
        'seqname'   => lc $seqname,
        'source'    => $source,
        'feature'   => $feature,
        'start'     => $start,
        'end'       => $end,
        'score'     => $score,
        'strand'    => $strand,
        'frame'     => $frame,
        'attribute' => $attribute
    };
}

__DATA__


__END__

=head1 NAME

 normalize_gff.pl - Normalize one GFF file by another

=head1 SYNOPSIS

 normalize_gff.pl [OPTION]... [[-i] FILE]...

 normalize_gff.pl -a exp.gff -b control.gff -c log2_ratio -f log2 -o norm_exp.gff

=head1 DESCRIPTION

 Long description

=head1 OPTIONS

 -a. gff-a
 -b, gff-b
 -c, scoring                    'log2_ratio' or 'difference'
 -f, feature
 -o, --output      <string>     output filename                          (STDOUT)
 -e, --error       <string>     output error filename                    (STDERR)
     --verbose     [integer]    print increasingly verbose error messages
     --quiet                    print no diagnostic or warning messages
     --version                  print current version
     --license                  print author's contact and copyright information
     --help                     print this information
     --manual                   print the plain old documentation page

=head1 VERSION

 0.0.1

=head1 REVISION

 $Rev: $:
 $Author: $:
 $Date: $:
 $HeadURL: $:
 $Id: $:

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
