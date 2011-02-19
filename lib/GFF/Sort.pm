package GFF::Sort;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;

use GFF qw/gff_to_string colname2num/;
use GFF::Parser;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(gff_sort);

sub gff_sort{
    my %opt = @_;
    my ($file,$overwrite,$colref,$manual) = @opt{'file','overwrite','cols','manual'};

    my @cols = @$colref;

    my $tmpfile = $file . ".sorted";
    croak "$file not read/writeable" unless (-r $file && -w $file);
    #croak "$tmpfile not writeable" unless (-r $tmpfile);

    my @sortexec = grep {-x} qw{/usr/bin/sort /bin/sort};
    if (!$manual && ($^O =~ /^(?:linux|darwin|freebsd)$/ && @sortexec)){
        my @colnums = map {colname2num $_} @cols;

        my $colstring = join ",", @colnums;

        my $sortbin = $sortexec[0];
        system("$sortbin -k $colstring -n $file > $tmpfile");
        rename $tmpfile, $file or croak ("couldn't rename $tmpfile to $file?") if $overwrite;
    }
    else {
        my $p = GFF::Parser->new(file => $file);
        my $arr = $p->slurp();

        my @sorted = sort {
            my $res;
            foreach my $col (@cols) {
                $res = $a->{$col} <=> $b->{$col};
                last if !$res;
            }
            $res;
        } @$arr;

        open my $fh, '>', $tmpfile;
        for (@sorted){
            say $fh gff_to_string($_);
        }
        close $fh;
        rename $tmpfile, $file or croak ("couldn't rename $tmpfile to $file?") if $overwrite;
    }
}

1;

