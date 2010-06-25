#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use version; our $VERSION = qv('0.0.1');
eval {
    require Term::ProgressBar::Simple;
};
croak "$0: requires Term::ProgressBar::Simple:$@%s" if $@;

GetOptions(
    \%ARGV,
    'input|i=s', 'output|o=s', 'error|e=s',
    'size|s=i', 'name|n=s',
    _meta_options( \%ARGV ),
) and $ARGV{input} or @ARGV or $ARGV{size} or pod2usage( -verbose => 1 );

my ( $INH, $OUTH, $ERRH ) = _prepare_io( \%ARGV, \@ARGV );

my ($name, $input_size) = (@ARGV{qw/name size/});

if (@ARGV) {
    $input_size =  0;
    $input_size += -s $_ for @ARGV;
}

use constant BUFFER_SIZE => 4096;
use constant UPDATE_PROB => 1/2;
my $buffer               =  q{};

my $progress = Term::ProgressBar::Simple->new( {
    name  => $name,
    count => $input_size,
    ETA   => 'linear',
    fh    => \$ERRH,
} );

if (@ARGV) {
    while (my $file = shift @ARGV) {
        open $INH, q{<}, $file or die qq{Cant read $file: $!};
        _cat( $INH, $OUTH );
    }
}
else {
    _cat( $INH, $OUTH );
}
close $INH;
close $OUTH;


sub _cat {
    my ($inh, $outh) = @_;
    while (my $read = sysread $inh, $buffer, BUFFER_SIZE) {
        if (syswrite $outh, $buffer) {
            $progress += $read unless int rand 1/UPDATE_PROB;
        }
    }
}

sub _meta_options {
    my ($opt) = @_;

    return (
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
    elsif (@$argv) { $INH = *STDIN }
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

 tpb.pl - Progress Viewer

=head1 SYNOPSIS

 tpb.pl [OPTION]... [[-i] FILE]...

=head1 DESCRIPTION

 Concatenate input to output with progress bar 

=head1 OPTIONS

 -n, --name        <string>     name of progress (useful when nesting multiple instances of tpb.pl)
 -s, --size        <integer>    size of input in bytes (ignored when not reading STDIN)
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
