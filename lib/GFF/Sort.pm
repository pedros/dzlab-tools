package GFF::Sort;
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use autodie;

use GFF::Parser;
use GFF::Misc;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(gff_is_sorted);
our @EXPORT = qw(gff_sort);

sub gff_sort{
    my %opt = @_;
    my ($file,$overwrite,$column,$manual) = @opt{'file','overwrite','column','manual'};
    
    croak 'need column!' unless $column;
    my $colnum = ($column eq 'start') ? 4 : ($column eq 'end') ? 5 : 
    croak ('gff_sort: column needs to be start or end');

    my $tmpfile = $file . ".sorted";
    croak "$file not read/writeable" unless (-r $file && -w $file);
    #croak "$tmpfile not writeable" unless (-r $tmpfile);

    my @sortexec = grep {-x} qw{/usr/bin/sort /bin/sort};
    if (!$manual && ($^O =~ /^(?:linux|darwin|freebsd)$/ && @sortexec)){
        
        my $sortbin = $sortexec[0];
        system("$sortbin -k $colnum -n $file > $tmpfile");
        rename $tmpfile, $file or croak ("couldn't rename $tmpfile to $file?") if $overwrite;
    }
    else {
        my $p = GFF::Parser->new(file => $file);
        my $arr = $p->slurp();

        my @sorted = sort {
            $a->get_column($column) <=> $b->get_column($column);
        } @$arr;

        open my $fh, '>', $tmpfile;
        for (@sorted){
            say $fh $_->to_string();
        }
        close $fh;
        rename $tmpfile, $file or croak ("couldn't rename $tmpfile to $file?") if $overwrite;
    }
    return $overwrite ? $file : $tmpfile;
}

sub gff_is_sorted{
    my ($file,$column) = @_;
    croak 'gff_is_sorted bad args' unless $file && $column;
    my $p = GFF::Parser->new(file => $file);
    my $first = $p->next();
    while (my $second = $p->next()){
        if ($first->$column >= $second->$column){
            return 0;
        }
        $first = $second;
    }
    return 1;
}

1;

