#!/usr/bin/env perl
# ___UNDOCUMENTED___

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use version; our $VERSION = qv('0.0.1');
use File::Basename;

use FindBin;
use lib "$FindBin::Bin/DZLab-Tools/lib";
use DZLab::Tools::RunUtils;

# use File::Spec;
#print map { File::Spec->rel2abs($_) ."\n" } @ARGV[0,1,2,3];exit;

GetOptions(
    \%ARGV,
    'input|i=s', 'output|o=s', 'error|e=s',
    _meta_options( \%ARGV ),
) and (@ARGV or $ARGV{input}) or pod2usage( -verbose => 1 );

my ( $INH, $OUTH, $ERRH ) = _prepare_io( \%ARGV, \@ARGV );

print $OUTH join ("\t", '# gene', map { $_ = basename $_; s/\..*//g; $_ } @{[@ARGV]}), "\n";

while (1) {

    my ($gene, @scores);

    my @lines = read_lines( @ARGV );
    last unless @lines;

    for ( @lines ) {
        my ($id, $sc, @fields);

        @fields = split /\t/, $_;
        @fields = split /\s+/, $_ unless @fields > 1;

        if (@fields == 2) {
            ($id, $sc) = @fields[0,1];
        }
        else {
            ($id) = $fields[-1] =~ m/ID=(\w+)/;
            ($sc) = $fields[5];
        }

        $id or die $fields[-1];
        defined $sc or die $_;

        $gene ||= $id;
        die "$gene different from $id in $_" unless $gene eq $id;
        $sc = 0 if q{.} eq $sc;
        push @scores, $sc;
    }

    print $OUTH join ("\t", $gene, @scores), "\n";
}

{
    my %handles;

    sub read_lines {
        my (@files, @lines) = @_;

        for (@files) {

            unless (exists $handles{$_}) {
                open my $IN, '<', $_ or die "Can't open $_ $!";
                $handles{$_} = $IN;
            }

            my $handle = $handles{$_};
            my $line = <$handle>;
            return unless defined $line;
            chomp $line;
            push @lines, $line;
        }
        return @lines;
    }
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
