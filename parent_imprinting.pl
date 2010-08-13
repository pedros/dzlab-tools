#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use version; our $VERSION = qv('0.0.1');
use List::Util qw/sum/;
use List::MoreUtils qw/uniq/;

my %statistics = (
    'CHI::phi'     => 1, 'CHI::tscore'   => 1, 'CHI::x2'           => 1,
    'Dice::dice'   => 1, 'Dice::jaccard' => 1,
    'Fisher::left' => 1, 'Fisher::right' => 1, 'Fisher::twotailed' => 1,
    'MI::ll'       => 1, 'MI::pmi'       => 1, 'MI::ps'            => 1, 'MI::tmi' => 1,
);

@ARGV{qw/statistics
         expected-maternal-imprinting-ratio
         expected-paternal-imprinting-ratio/}
= ( 'Fisher::twotailed', 2/3, 1/3 );

@ARGV{qw/maternal-test paternal-test
         maternal-control paternal-control/}
= ([], [], [], []);

GetOptions(
    \%ARGV,
    'input|i=s', 'output|o=s', 'error|e=s',
    'maternal-test|mt=i{,}', 'maternal-control|mc=i{,}',
    'paternal-test|pt=i{,}', 'paternal-control|pc=i{,}',
    'expected-maternal-imprinting-ratio|m=f',
    'expected-paternal-imprinting-ratio|p=f',
    'statistics|s=s',
    _meta_options( \%ARGV ),
) and (@ARGV or $ARGV{input}) or pod2usage( -verbose => 1 );

my ( $INH, $OUTH, $ERRH ) = _prepare_io( \%ARGV, \@ARGV );

if (exists $statistics{$ARGV{statistics}}) {
    eval "use Text::NSP::Measures::2D::$ARGV{statistics}";
}
else {
    pod2usage( -sections => ['OPTIONS'],
               -message  => "Invalid statistic: $ARGV{statistics}",
               -verbose  => 99 );
}

chomp( my @headers = split /\t/, scalar <$INH> );

my @uniq_sorted_headers
= sort { $a <=> $b }
uniq
grep { defined $headers[$_] }
map { ref $_ ? @$_ : $_ }
@ARGV{qw/maternal-test paternal-test maternal-control paternal-control/};

print $OUTH join( "\t",
            $headers[0],
            @headers[@uniq_sorted_headers],
            $ARGV{statistics},
            $ARGV{debug} ? qw/mt mc pt pc/ : (),
        ), "\n";

while ( <$INH> ) {
    # columns: mm, mp, pp. pm
    # contigency table:
    #  |test   |control|
    # -|-------|-------|----
    # m| n11   | n12   | n1p
    # -|-------|-------|----
    # p| n21   | n22   |
    #  |-------|-------|----
    #  | np1   |       | npp 
    chomp;
    my @fields = split /\t/;

    my ($n11, $n12, $n21, $n22);

    $n11 = sum @fields[ @{ $ARGV{'maternal-test'} } ];
    $n12 = sum @fields[ @{ $ARGV{'paternal-test'} } ];

    if (grep {@$_} @ARGV{qw/maternal-control paternal-control/}) {
        $n21 = sum @fields[ @{ $ARGV{'maternal-control'} } ];
        $n22 = sum @fields[ @{ $ARGV{'paternal-control'} } ];
    }
    else {
        $n21 = $ARGV{'expected-maternal-imprinting-ratio'} * ($n11 + $n12);
        $n22 = $ARGV{'expected-paternal-imprinting-ratio'} * ($n11 + $n12);
    }

    die Dumper $n11, $n12, \%ARGV unless defined $n21 and defined $n22;

    my $p = get_p( (map { sprintf "%.0f", $_ } $n11, $n12, $n21, $n22), $fields[0]);

    print $OUTH join( "\t",
                $fields[0],
                @fields[@uniq_sorted_headers],
                sprintf( "%g", $p),
                $ARGV{debug} ? map {sprintf "%g",  $_}( $n11, $n12, $n21, $n22 ) : (),
            ), "\n";
}

sub get_p {
    my ($n11, $n12, $n21, $n22, $gene) = @_;

    if ($n11 + $n12 + $n21 + $n22) {
        return calculateStatistic
        (
            n11 => $n11,
            n1p => ($n11 + $n12),
            np1 => ($n11 + $n21),
            npp => ($n11 + $n12 + $n21 + $n22)
        );
    }

   if ( my $error_code = getErrorCode() ) {
        print $ERRH "$error_code: ", getErrorMessage(), "\n";
    }
    else {return 0}
}

sub _meta_options {
    my ($opt) = @_;

    return (
        'debug',
        'quiet'     => sub { $opt->{quiet}   = 1;          $opt->{verbose} = 0 },
        'verbose:i' => sub { $opt->{verbose} = $_[1] // 1; $opt->{quiet}   = 0 },
        'version'   => sub { pod2usage( -sections => ['VERSION', 'REVISION'],
                                        -verbose  => 99 )                      },
        'license'   => sub { pod2usage( -sections => ['AUTHOR', 'COPYRIGHT'],
                                        -verbose  => 99 )                      },
        'usage'     => sub { pod2usage( -sections => ['SYNOPSIS'],
                                        -verbose  => 99 )                      },
        'help'      => sub { pod2usage( -verbose  => 1  )                      },
        'manual'    => sub { pod2usage( -verbose  => 2  )                      },
    );
}

sub _prepare_io {
    my ($opt, $argv) = @_;

    my ($INH, $OUTH, $ERRH);
    
    # If user explicitly sets -i, put the argument in @$argv
    unshift @$argv, $opt->{input} if exists $opt->{input};

    # Allow in-situ arguments (equal input and output filenames)
    if (    exists $opt->{input} and exists $opt->{output}
               and $opt->{input} eq $opt->{output} ) {
        open $INH, q{<}, $opt->{input}
            or croak "Can't read $opt->{input}: $!";
        unlink $opt->{output};
    }
    else { $INH = *ARGV }

    # Redirect STDOUT to a file if so specified
    if ( exists $opt->{output} and q{-} ne $opt->{output} ) {
        open $OUTH, q{>}, $opt->{output}
            or croak "Can't write $opt->{output}: $!";
    }
    else { $OUTH = *STDOUT }

    # Log STDERR if so specified
    if ( exists $opt->{error} and q{-} ne $opt->{error} ) {
        open $ERRH, q{>}, $opt->{error}
            or croak "Can't write $opt->{error}: $!";
    }
    elsif ( exists $opt->{quiet} and $opt->{quiet} ) {
        use File::Spec;
        open $ERRH, q{>}, File::Spec->devnull
            or croak "Can't write $opt->{error}: $!";
    }
    else { $ERRH = *STDERR }

    return ( $INH, $OUTH, *STDERR = $ERRH );
}

__DATA__


__END__

=head1 NAME

 APerlyName.pl - Short description

=head1 SYNOPSIS

 APerlyName.pl [OPTION]... [[-i] FILE]...

=head1 DESCRIPTION

 Long description

=head1 OPTIONS

 -i, --input       <string>     input filename                           (STDIN)
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
