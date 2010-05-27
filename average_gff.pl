#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;
use version; our $VERSION = qv('0.0.1');

use List::Util qw(sum min);

GetOptions(
    \%ARGV,
    'input|i=s', 'output|o=s', 'error|e=s',
    'min-meth|m=f', 'max-meth|M=f', 'loci-list|l=s', 'loci-tag|t=s',
    _meta_options( \%ARGV ),
) and (@ARGV or $ARGV{input}) or pod2usage( -verbose => 1 );

my ( $INH, $OUTH, $ERRH ) = _prepare_io( \%ARGV, \@ARGV );

my @totals;

$ARGV{'loci-list'} = index_list( @ARGV{qw/loci-list loci-tag/} )
if exists $ARGV{'loci-list'};

GFF:
while ( <$INH> ) {
    chomp;
    my ($score, $attribute) = (split /\t/, $_)[5, 8];

    exists $ARGV{'min-meth'}  and $score < $ARGV{'min-meth'} and next GFF;
    exists $ARGV{'max-meth'}  and $score > $ARGV{'max-meth'} and next GFF;

    my %attributes = $attribute =~ m/(\w+)=([^;]+)/g;

    exists $ARGV{'loci-list'} and not exists $ARGV{'loci-list'}->{$attributes{ID}}
    and next GFF;

    delete $attributes{ID};

    $attributes{score} = $score;

    push @totals, \%attributes;
}

my $totals = squash( \@totals );

print $OUTH join ("\n",
                  "##methylation:\t" . sprintf ("%.5g", $totals->{score}),
                  "##n.loci:\t$totals->{sites}",
                  "##sites:\t$totals->{n}"
              ), "\n";

delete @{$totals}{qw/score sites n/};

print $OUTH join ("\t", '#pos', qw/a c g t/), "\n";

for my $position (sort {$a <=> $b} keys %$totals) {

    if ($position =~ m/\d+/) {
        print $OUTH $position,
                    "\t",
                    join ("\t",
                          map { sprintf "%.5g", $_ } @{$totals->{$position}}
                    ),
                    "\n";
    }
}

sub index_list {
    my ($list, $tag) = @_;

    $tag ||= '';
    my $counter = 0;

    my %list = ();
    open my $LIST, '<', $list or croak "Can't open $list for reading";
    while (<$LIST>) {

        my ($id, $freq, $alt) = (0, 0, 0);
        my @fields = split /\t/, $_;

        if (@fields < 9) {
            ($id, $freq, $alt) = @fields;
        }
        elsif (@fields == 9) {
            ($id, $freq, $alt) = @fields[8, 5, 0];
            $id =~ s/\W*$tag=?(\w+).*/$1/;
        }
        else {
            ($id) = $fields[0];
        }
        $id =~ s/[\r\n]//g;
        $list{$id} = [$freq, $alt, $counter++];
    }
    close $LIST or carp "Can't close $list after reading";

    return \%list;
}

sub squash {
    my ($totals_ref) = @_;

    my %totals;

    my $total_n  = 0;
       $total_n += sum $_->{n} for @$totals_ref;

    for my $loci (@$totals_ref) {   

        my $weight = $loci->{n} / $total_n;

        while (my ($k, $v) = each %$loci) {

            my @v = split /,/, $v;

            if ('n' eq $k) {
                $totals{n} += $v;
            }
            else {
     
                if (4 < @v) {
                    my $min = min @v;
                    $min = do { $_ if $v[$_] eq $min} for 0 .. @v - 1;
                    splice @v, $min, 1;
                }

                $totals{$k}[$_] += ($v[$_] * $weight)
                for 0 .. @v - 1;
            }

        }
    }

    $totals{sites} = @$totals_ref;
    $totals{score} = $totals{score}->[0];
    return \%totals;
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
