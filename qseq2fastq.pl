#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use version; our $VERSION = qv('0.0.1');

GetOptions( 
    \%ARGV, 
    'input|i=s{,}', 'output|o=s', 'error|e=s',
    _meta_options( \%ARGV ),
)
and ( @ARGV or $ARGV{input} )
or pod2usage( -verbose => 1 );

my ( $INH, $OUTH, $ERRH ) = _prepare_io( \%ARGV, \@ARGV );

use constant {
    MACHINE_ID     => 0,
    RUN_NUMBER     => 1,
    LANE_NUMBER    => 2,
    TILE_NUMBER    => 3,
    X_COORDINATE   => 4,
    Y_COORDINATE   => 5,
    INDEX          => 6,
    READ           => 7,
    SEQUENCE       => 8,
    QUALITY_SCORES => 9,
    FLAG_FILTER    => 10,
};

while (<$INH>) {
    chomp;

    print $OUTH fastq( split /\t/, $_ );
}


sub fastq {
    my (@read) = @_;

    $read[SEQUENCE] =~ s{\.}{N}g;

    my $fastq = join(
        q{:},
        join( q{_}, @read[MACHINE_ID,RUN_NUMBER] ),
        @read
        [
            LANE_NUMBER,
            TILE_NUMBER,
            X_COORDINATE,
            Y_COORDINATE,
        ],
    );

    return
    sprintf "\@%s\n%s\n+%s\n%s\n",
    $fastq, $read[SEQUENCE], $fastq, $read[QUALITY_SCORES];
}

sub _meta_options {
    my ($opt) = @_;

    return (
        'quiet' => sub { $opt->{quiet} = 1; $opt->{verbose} = 0 },
        'verbose:i' => sub { $opt->{verbose} = $_[1] // 1; $opt->{quiet} = 0 },
        'version' => sub { pod2usage( -sections => [ 'VERSION', 'REVISION' ],
                -verbose  => 99)
        },
        'license' => sub {
            pod2usage(
                -sections => [ 'AUTHOR', 'COPYRIGHT' ],
                -verbose  => 99
            )
        },
        'usage' => sub {
            pod2usage(
                -sections => ['SYNOPSIS'],
                -verbose  => 99
            )
        },
        'help'   => sub { pod2usage( -verbose => 1 ) },
        'manual' => sub { pod2usage( -verbose => 2 ) },
    );
}

sub _prepare_io {
    my ( $opt, $argv ) = @_;

    my ( $INH, $OUTH, $ERRH );

    # If user explicitly sets -i, put the argument in @$argv
    unshift @$argv, $opt->{input} if exists $opt->{input};

    # Allow in-situ arguments (equal input and output filenames)
    if (    exists $opt->{input}
        and exists $opt->{output}
        and $opt->{input} eq $opt->{output} )
    {
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


=head1 NAME

 qseq2fastq.pl - Convert qseq to fastq

=head1 SYNOPSIS

 qseq2fastq.pl [OPTION]... [[-i] FILE]...

 qseq2fastq.pl s_1_1_*.qseq > s_1_1_sequence.txt
 qseq2fastq.pl -i s_2_*.qseq -o s_2_sequence.fastq
 qseq2fastq.pl --input s_3_*.qseq --output s_3_sequence.fastq


=head1 DESCRIPTION

 Long description

=head1 OPTIONS

 -i, --input       <string>     input filename(s)                        (STDIN)
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
