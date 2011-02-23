package GFF::Split;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;
use File::Basename;
use GFF::Parser;
use Cwd;
use File::Spec;
use File::Temp qw/tempdir/;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(gff_split);

# return hash of sequence/feature to filename
sub gff_split{
    my (%opt) = @_;
    my ($file, $sequence,$feature) = @opt{'file', 'sequence','feature'};

    croak "gff_split(file => 'filename',[sequence|feature] => 'all',tmpdir => 1')" unless ($sequence xor $feature);

    my ($name, $path, $suffix);
    if (ref $file eq 'GLOB'){
        $name = 'FILEHANDLE';
        $path = getcwd;
        $suffix = '.gff';
    } else{
        ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);
    }

    my $outdir = $opt{tmpdir} ? tempdir(CLEANUP => 1) : $path; # if directory given, use it-- otherwise use same dir

    my %file_handle;
    my %files;
    my $parser = GFF::Parser->new(file => $file);
    RECORD:
    while (my $gff = $parser->next()){
        my $current = $sequence ? $gff->sequence : $gff->feature;
        next RECORD unless defined $current;

        if (($feature && ($feature eq 'all' || $feature eq $gff->feature))
            ||
            ($sequence && ($sequence eq 'all' || $sequence eq $gff->sequence)))
        {
            my $out = File::Spec->catfile($outdir, $name) . "-$current" . $suffix;

            if (! exists $file_handle{$current}){
                open $file_handle{$current}, '>', $out
                    or croak "Can't open $out for writing: $!";
                $files{$current} = $out;
            }

            say {$file_handle{$current}} $gff->to_string();
        }
    }
    foreach my $fh (values %file_handle) {
        close $fh;
    }
    return \%files;
}

1;

