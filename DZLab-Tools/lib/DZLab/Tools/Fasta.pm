package DZLab::Tools::Fasta;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(slurp_fasta);


=head1 slurp_fasta
return hash of seqid (the first word after the '>') to the sequence.
=cut
sub slurp_fasta {
    my ($file, $opts) = @_;
    return {} unless $file;

    $opts->{-l} ||= 0;

    my %accum = ();

    open my $fh, '<', $file or croak "Can't open $file: $?";

    my $current;
    my @buffer;
    while (defined(my $line = <$fh>)) {
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

    close $fh;

    return \%accum;
}
1;


