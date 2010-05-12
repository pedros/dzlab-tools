#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;

my $INH  = *ARGV;
my $ERRH = *STDERR;
my $OUTH = *STDOUT;
GetOptions (
    \%ARGV,
    'input|i=s',
    'output|o=s',
    'debug:i',
    'quiet'           => sub {$ARGV{quiet}   = 1; no  diagnostics; no warnings },
    'verbose'         => sub {$ARGV{verbose} = 1; use diagnostics; use warnings},
    'version'         => sub {pod2usage -sections => ['VERSION','REVISION'],
                                        -verbose  => 99                       },
    'license'         => sub {pod2usage -sections => ['AUTHOR','COPYRIGHT'],
                                        -verbose  => 99                       },
    'usage'           => sub {pod2usage -sections => ['SYNOPSIS'],
                                        -verbose  => 99                       },
    'help'            => sub {pod2usage -verbose  => 1                        },
    'manual'          => sub {pod2usage -verbose  => 2                        },
    )
    or pod2usage (-verbose => 1);


 IO:
{
    # We use the ARGV magical handle to read input
    # If user explicitly sets -i, put the argument in @ARGV
    if (exists $ARGV{input}) {
        unshift @ARGV, $ARGV{input};
    }

    # Allow in-situ arguments (equal input and output filenames)
    # FIXME: infinite loop. Why?
    if (exists $ARGV{input} and exists $ARGV{output}
        and $ARGV{input} eq $ARGV{output}) {
        croak "Bug: don't use in-situ editing (same input and output files";
        open $INH, q{<}, $ARGV{input}
        or croak "Can't read $ARGV{input}: $!";
        unlink $ARGV{input};
    }

    # Redirect STDOUT to a file if so specified
    if (exists $ARGV{output} and q{-} ne $ARGV{output}) {
        open $OUTH, q{>}, $ARGV{output}
        or croak "Can't write $ARGV{output}: $!";
    }
}

my @buffer;
my %chr_counts;

MAIN:
while(<$INH>){
    chomp;
    local $_ = [split /\t/, $_];

    # flush buffer
    if (@buffer and ("@{$_}[0,3,4,6]" ne "@{$buffer[-1]}[0,3,4,6]")) {
        $buffer[-1]->[5] = scalar @buffer - 1;
        $chr_counts{$buffer[-1]->[0]} += $buffer[-1]->[5];
        print {$OUTH} join ("\t", @{$buffer[-1]}), "\n";
        @buffer = ();
    }
    
    push @buffer, $_;
}

print {$OUTH} join ("\t", @{$buffer[-1]}), "\n"
if @buffer;

for (sort keys %chr_counts) {
    print {$ERRH} "$_\t$chr_counts{$_}\n";
}

__DATA__

__END__
=head1 NAME

 filter_repeats.pl - Filter out GFF records with same seqid, start, end, strand

=head1 SYNOPSIS

 filter_repeats.pl [OPTION]... [FILE]...

 perl filter_repeats.pl some.gff -o filtered.gff 2> some.repeats

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
