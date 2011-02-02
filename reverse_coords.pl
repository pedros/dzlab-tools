
# ___UNDOCUMENTED___
#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Data::Dumper;

my $gfffile = $ARGV[0];
my $referencefile = $ARGV[1];

# holds name of chromosomes as keys and length of chromosomes in bp as values
my %reference = ();

# begins limited scope for extracting names and lengths of chromosomes in reference file
{
    # reads in the reference genome file into @fastaseq
    open my $REF, '<', "$referencefile" or croak "Can't open file: $referencefile";
    my @fastaseq = <$REF>;
    close $REF;

    # find and store indices for each chromosome change and corresponding descriptions
    my ( @idx, @dsc ) = ();
    for my $i ( 0 .. @fastaseq - 1 ) { ### Indexing $referencefile...  % done
        if ( $fastaseq[$i] =~ m/^>/ ) {
            $fastaseq[$i] =~ s/>//g;
            $fastaseq[$i] = ( split /\s/, "$fastaseq[$i]" )[0];
            $fastaseq[$i] =~ tr/A-Z/a-z/;
            push @idx, $i;
            push @dsc, $fastaseq[$i];
        }
    }

    # gets and saves each chromosome's sequence and reverse complemented sequence
    for my $j ( 0 .. @idx - 1 ) { ### Loading $referencefile into memory...  % done
        my $line;
        if ( $j == scalar @idx - 1 ) {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. @fastaseq - 1]);
        }
        else {
            $line = join( q{}, @fastaseq[ $idx[$j] + 1 .. $idx[$j + 1] - 1]);
        }
        $line =~ s/[\n\r]//g;
        $reference{ $dsc[$j] } = length $line;
    }
}

open my $GFF, '<', $gfffile or croak "can't open $gfffile";

GFF_READ_LOOP:
while (my $gff_line = <$GFF>) {
    next GFF_READ_LOOP if ($gff_line =~ m/^#.*$|^\s*$/);

    my @fields =  split /\t/, $gff_line;
    $fields[0] =~ tr/A-Z/a-z/;
    $fields[3] =  $reference{$fields[0]} + 1 - $fields[3];
    $fields[4] =  $reference{$fields[0]} + 1 - $fields[4];

    ##temporary
    $fields[6] = q{-};

    print join "\t", @fields;
}

close $GFF;
