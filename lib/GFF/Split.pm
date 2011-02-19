package GFF::Split;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use File::Basename;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(gff_split);

use FindBin;
use lib "$FindBin::Bin/../lib";
use GFF::Parse qw/gff_to_string gff_make_iterator/;

# return hash of sequence/feature to filename
sub gff_split{
    my ($file, %opt) = @_;
    my ($seqname,$feature) = @opt{'seqname','feature'};

    croak "gff_split('filename',[sequence|feature] => 'all',directory => '/tmp')" unless ($seqname xor $feature);

    my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);

    my $outdir = $opt{directory} // $path; # if directory given, use it-- otherwise use same dir

    my %file_handle;
    my %files;
    my $iter = gff_make_iterator(file => $file);
    RECORD:
    while (defined(my $gff = $iter->())){
        next unless ref $gff eq 'HASH';

        my $current = $seqname ? $gff->{seqname} : $gff->{feature};
        next RECORD unless defined $current;

        if (($feature && ($feature eq 'all' || $feature eq $gff->{feature}))
            ||
            ($seqname && ($seqname eq 'all' || $seqname eq $gff->{seqname})))
        {
            my $out = $outdir . $name . "-$current" . $suffix;

            if (! exists $file_handle{$current}){
                open $file_handle{$current}, '>', $out
                    or croak "Can't open $out for writing: $!";
                $files{$current} = $out;
            }

            say {$file_handle{$current}} gff_to_string $gff;
        }
    }
    foreach my $fh (values %file_handle) {

    }
    return \%files;
}

1;

__END__

=head1 NAME

 split_gff - Split GFF files by sequence ID or feature

=head1 SYNOPSIS

 split_gff.pl -f exon all_features.gff  # filter by exon, create file all_features-exon.gff
 split_gff.pl -s chr1 all_sequences.gff # filter by chromosome, create file all_sequences-chr1.gff
 split_gff.pl -f all all_features.gff   # create multiple files, one per feature

=head1 DESCRIPTION

=head1 OPTIONS

 split_gff.pl [OPTION]... [FILE]...

 -f, --feature     feature used to filter GFF file by ('all' generates one file per feature)
 -s, --sequence    sequence ID used to filter GFF file by ('all' generates one file per sequence ID)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: 249 $:
 $Author: psilva $:
 $Date: 2010-01-11 21:24:34 -0800 (Mon, 11 Jan 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/split_gff.pl $:
 $Id: split_gff.pl 249 2010-01-12 05:24:34Z psilva $:

=head1 AUTHOR

 Pedro Silva <psilva@nature.berkeley.edu/>
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
