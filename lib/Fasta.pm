package Fasta;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(slurp_fasta format_fasta count_fasta count_fasta_complete);

=head1 EXPORTED FUNCTIONS

=head2 slurp_fasta

return hash of seqid (the first word after the '>') to the sequence.

=cut
sub slurp_fasta {
    my ($file, $opts) = @_;
    return {} unless $file;

    $opts->{-l} //= 0;

    my %accum = ();

    open my $fh, '<', $file or croak "Can't open $file: $?";

    my $current;
    my @buffer;
    while (defined(my $line = <$fh>)) {
        $line =~ tr/\r\n//d;
        if ($line =~ /^>(\w+)/){
            if ($current){
                $accum{$current} = join q{}, @buffer;
                @buffer = ();
            }
            $current = lc $1;
        } 
        else{
            chomp $line;
            push @buffer, $opts->{-l} ? lc $line : $line;
        }
    }

    $accum{$current} = join q{}, @buffer if ($current) ; # last seq
    $accum{$current} =~ s/\s*//g;

    close $fh;

    return \%accum;
}

sub count_fasta_complete {
    my ($file) = @_;
    return {} unless $file;

    my %accum = ();

    open my $fh, '<', $file or croak "Can't open $file: $?";

    my $current;
    while (defined(my $line = <$fh>)) {
        $line =~ tr/\r\n//d;
        if ($line =~ /^>(\w+)/){
            $current = lc $1;
        } 
        else{
            my ($len, $bp) = (length $line, $line =~ tr/acgtACGT//);
            $accum{$current}{length} += $len;
            $accum{$current}{bp}     += $bp;
            $accum{$current}{nonbp}  += $len-$bp;
        }
    }
    close $fh;

    return \%accum;
}
sub count_fasta {
    my ($file) = @_;
    my $counts = count_fasta_complete($file);
    return {map { $_ => $counts->{$_}{length} } keys %$counts};
}


=head2 format_fasta

format a header and seq for printing as a fasta.

=cut 
sub format_fasta{
    my ($header, $seq, $width) = @_;

    croak "need a sequence and header" unless ($seq and $header);

    my @buffer;
    $buffer[0] = ">$header" if defined $header;
    $width //= 80;
    my $length = length $seq;

    for (my $position = 0; $position < $length; $position += $width){
        push @buffer, substr $seq, $position, $width;
    }
    return (join "\n", @buffer) . "\n";
}

1;


