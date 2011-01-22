#!/usr/bin/env perl -w
# parse_ecker.pl --- Parses GSM276809.txt data set from Ecker et all, Cell 133
# Author: Pedro Silva <psilva@dzlab.pmb.berkeley.edu>
# Created: 20 Jan 2009
# Version: 0.01

use warnings;
use strict;
use Carp;
use Smart::Comments;

# begins limited scope for extracting names and lengths of chromosomes in reference file
sub index_fasta {
    my $referencefile = shift;

    my %reference = ();
    
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
            $line = join( q{}, 'NN', @fastaseq[ $idx[$j] + 1 .. @fastaseq - 1] , 'NN');
        }
        else {
            $line = join( q{}, 'NN', @fastaseq[ $idx[$j] + 1 .. $idx[$j + 1] - 1] , 'NN');
        }
        $line =~ s/[\n\r]//g;
        $reference{ $dsc[$j] }     = length $line;
        $reference{"$dsc[$j]-seq"} = $line;
#        $reference{"$dsc[$j]-rc"}  = reverse_complement ($line);
    }
    return \%reference;
}


# Converts input sequence to reverse complement
# Returns scalar $string with processed sequence
sub reverse_complement {
    my $shortseq = shift;
    $shortseq =~ tr/ACGTacgt/TGCAtgca/;
    $shortseq =~ s/\n//;
    return scalar reverse $shortseq;
}


## MAIN PROGRAM

my %genome_ref = %{index_fasta ($ARGV[0])};

open my $EXT, '<', $ARGV[1] or croak "Can't open $ARGV[1]";

INPUT:
while (<$EXT>) { ### Parsing $ARGV[1]...
    $_ =~ s/[\r\n]//g;
    next INPUT if $. == 1;
    my @fields  = split (/\t/, $_) if $. > 1;

    next INPUT if ($fields[0] !~ m/[CM]/i or $fields[1] eq q{-});
    
    my $read_len = length $fields[4];
    $fields[0] =~ tr/A-Z/a-z/;
    
    if ($fields[1] eq q{+}) {

        print join("\t",
                   "Chr$fields[0]",
                   '.',
                   '/1:' . $fields[4],
                   $fields[2],
                   $fields[3],
                   1,
                   $fields[1],
                   '.',
                   'target=' . substr ($genome_ref{"chr$fields[0]-seq"}, $fields[2] - 1, $read_len + 4),
                );                       
    }
    else {

        my $reverse_seq = reverse $fields[4];
            
        print join("\t",
                    "Chr$fields[0]",
                    '.',
                    '/1:' . $reverse_seq,
                    $fields[2],
                    $fields[3],
                    1,
                    $fields[1],
                    '.',
                    'target=' . reverse_complement (substr ($genome_ref{"chr$fields[0]-seq"}, $fields[2] - 1, $read_len + 4))
                );
    }
    print "\n";
}

close $EXT;


__END__

=head1 NAME

parse_ecker.pl - Describe the usage of script briefly

=head1 SYNOPSIS

parse_ecker.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for parse_ecker.pl, 

=head1 AUTHOR

Pedro Silva, E<lt>psilva@dzlab.pmb.berkeley.eduE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 by Pedro Silva

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
