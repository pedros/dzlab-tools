#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config gnu_getopt no_bundling);
use Pod::Usage;

my %statistics = (
    'CHI::phi'          => 1,
    'CHI::tscore'       => 1,
    'CHI::x2'           => 1,
    'Dice::dice'        => 1,
    'Dice::jaccard'     => 1,
    'Fisher::left'      => 1,
    'Fisher::right'     => 1,
    'Fisher::twotailed' => 1,
    'MI::ll'            => 1,
    'MI::pmi'           => 1,
    'MI::ps'            => 1,
    'MI::tmi'           => 1,
);

GetOptions (
    \%ARGV,
    'maternal-proportion|m=f',
    'columns|c=i{4}' => \@{$ARGV{columns}}, # mm, mp, pm, pp
    'statistics|s=s' => sub {exists $statistics{$_[1]}
                             ? do {
                                 $ARGV{$_[0]} = $_[1];
                                 eval "use Text::NSP::Measures::2D::$_[1]";
                             }
                             : pod2usage -sections => ['OPTIONS'],
                                         -message  => "Invalid statistic: $_[1]",
                                         -verbose  => 99                       },
    # meta-options
    'input|i=s',
    'output|o=s',
    'debug:i',
    'quiet'          => sub {$ARGV{quiet}   = 1; no  diagnostics; no warnings  },
    'verbose'        => sub {$ARGV{verbose} = 1; use diagnostics; use warnings},
    'version'        => sub {pod2usage -sections => ['VERSION','REVISION'],
                                        -verbose  => 99                        },
    'license'        => sub {pod2usage -sections => ['AUTHOR','COPYRIGHT'],
                                        -verbose  => 99                        },
    'usage'          => sub {pod2usage -sections => ['SYNOPSIS'],
                                        -verbose  => 99                        },
    'help'           => sub {pod2usage -verbose  => 1                          },
    'manual'         => sub {pod2usage -verbose  => 2                          },
    )
    or pod2usage (-verbose => 1);

my $INH  = *ARGV;
my $ERRH = *STDERR;
my $OUTH = *STDOUT;

# We use the ARGV magical handle to read input
# If user explicitly sets -i, put the argument in @ARGV
if (exists $ARGV{input}) {
    unshift @ARGV, $ARGV{input};
}
# Allow in-situ arguments (equal input and output filenames)
# FIXME: infinite loop. Makes sense, but why does Conway recommend it?
if (exists $ARGV{input} and exists $ARGV{output}
    and    $ARGV{input} eq         $ARGV{output}) {
    croak "Bug: don't use in-situ editing (same input and output files"; #
    open $INH, q{<}, $ARGV{input} or croak "Can't read $ARGV{input}: $!";
    unlink $ARGV{input};
}
# Redirect STDOUT to a file if so specified
if (exists $ARGV{output} and q{-} ne $ARGV{output}) {
    open $OUTH, q{>}, $ARGV{output}
    or croak "Can't write $ARGV{output}: $!";
}

MAIN:
my @headers = split /\t/, scalar <$INH>;
chomp $headers[-1];

print join("\t",
           @headers[0,1,2], $ARGV{statistics},
           @headers[3,4],   $ARGV{statistics}
       ), "\n";

while(<$INH>){
    chomp;
    my @values = split /\t/, $_;

    # columns: mm, mp, pp. pm
    # contigency table:
    #  |observ|expect|
    # -|------|------|----
    # m| n11  | n12  | n1p
    # -|------|------|----
    # p| n21  | n22  |
    #  |------|------|----
    #  | np1  |      | npp 

    my $n11 =  $values[$ARGV{columns}->[0]];
    my $n21 =  $values[$ARGV{columns}->[1]];

    my $n12 = ($n11 + $n21)
              * $ARGV{'maternal-proportion'};
    my $n22 = ($n11 + $n21)
              * (1 - $ARGV{'maternal-proportion'});

    print join("\t",
               @values[0,1,2],
               sprintf("%g", get_p($n11, $n12, $n21, $n22))
           ), "\t";

    $n11 =  $values[$ARGV{columns}->[2]];
    $n21 =  $values[$ARGV{columns}->[3]];

    $n12 = ($n11 + $n21)
              * (1 - $ARGV{'maternal-proportion'});
    $n22 = ($n11 + $n21)
              * $ARGV{'maternal-proportion'};    

    print join("\t",
               @values[3,4],
               sprintf("%g", get_p($n11, $n12, $n21, $n22))
           ), "\n";
}

## done
sub get_p {
    my ($n11, $n12, $n21, $n22) = @_;

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
        print STDERR "$error_code: ", getErrorMessage(), "\n";
    }
    else {return 0}
}


__DATA__

__END__

=head1 NAME

 parent_imprinting.pl - Short description

=head1 SYNOPSIS

 parent_imprinting.pl [OPTION]... [FILE]...

=head1 DESCRIPTION

 Long description

=head1 OPTIONS

 -i, --input       filename to read from                            (STDIN)
 -o, --output      filename to write to                             (STDOUT)
     --debug       print additional information
     --verbose     print diagnostic and warning messages
     --quiet       print no diagnostic or warning messages
     --version     print current version
     --license     print author's contact and copyright information
     --help        print this information
     --manual      print the plain old documentation page

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
