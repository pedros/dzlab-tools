#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use feature 'switch';
use autodie;
use Pod::Usage;
use Getopt::Long;

my $help;
my $input;
my $output;
my $rep;

my $result = GetOptions (
    "input|i=s"       => \$input,
    "output|o=s"      => \$output,
    "replacement|r=s" => \$rep,
    "help"            => \$help,
);
pod2usage(-verbose => 99, -sections => [qw/NAME SYNOPSIS OPTIONS/]) if (!$result || !$input || !$rep);

my $ofh;
given ($output){
    when (q{-}){
        $ofh = \*STDOUT;
    }
    when (undef){
        my $orig = $input;
        $input .= '.bak';
        rename $orig, $input;
        open $ofh, '>', $orig;
    }
    default {
        open $ofh, '>', $output;
    }
}

open my $ifh, q{<}, $input;

while (defined(my $line = <$ifh>)){
    my @split = split /\t/, $line;
    die "line $. does not look like a GFF line." 
    unless @split == 9;
    $split[2] = $rep;
    print $ofh join qq{\t}, @split;
}

close $ifh;
close $ofh;

=head1 NAME

rename_feature.pl - replace the feature column (column 3) 

=head1 SYNOPSIS

rename the windows in genes.gff to "new_window_name", creating a backup in genes.gff.bak

 rename_feature.pl -i genes.gff -r new_window_name 

create new file genes-new.gff with new_window_name as feature.

 rename_feature.pl -i genes.gff -r new_window_name -o genes-new.gff

=head1 OPTIONS

=over

=item --input <file> | -i <file>

Input file.  

=item --output <file> | -o <file>

Output file.  If this is omitted, input file will be overwritten AND a backup will be saved with a '.bak' file
extension.  If '-', output to screen.

=item --replacement "new name" | -r "new name"

New name for feature column.

=back

=cut

