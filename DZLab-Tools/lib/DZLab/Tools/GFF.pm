package DZLab::Tools::GFF;

use strict; use warnings;
use version; our $VERSION = '0.0.1';

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(gff_read gff_make_iterator);

use Carp;


sub gff_read {
    my ($gff_line) = @_;

    # get pragmas of the form: '##gff-version 3' and '##sequence-region ctg123 1 1497228'
    return [split /\s+/, $1] if $gff_line =~ m/^
                                               \s*
                                               \#{2}
                                               \s*
                                               (.*)
                                               $/mx;

    # ignore blank lines and lines starting with '#'
    return [] if $gff_line =~ m/^ \s* (?:\#+ .*)? $/mx;

    my ($seqname, $source, $feature, $start, $end,
        $score,   $strand, $frame,   $attribute
    ) = split m/\t/xm, $gff_line || return;

    $attribute =~ s/[\r\n]//mxg;

    my %attributes = map { split /=/, $_ } split /;/, $attribute;
    @attributes{keys %attributes} = map { /,/ ? [split /,/] : $_ } values %attributes;

    return {
        seqname   => lc $seqname,
        source    => $source,
        feature   => $feature,
        start     => $start,
        end       => $end,
        score     => $score,
        strand    => $strand,
        frame     => $frame,
        attribute => $attribute,
        attributes=> \%attributes
    };
}

sub gff_make_iterator {
    my %options = @_;

    my $parser     = $options{parser};
    my $file       = $options{file};
    my $GFF_HANDLE = $options{handle};

    croak
        "Need parser function reference and file name or handle to build iterator"
        unless $parser
            and ref $parser eq 'CODE'
            and (
                ( defined $file and -e $file )
                xor(defined $GFF_HANDLE and ref $GFF_HANDLE eq 'GLOB'
                    or $GFF_HANDLE eq 'ARGV'
                )
            );

    if ($file) {
        open $GFF_HANDLE, '<', $file
            or croak "Can't read $file: $!";
    }

    return sub {
        $parser->( scalar <$GFF_HANDLE> );
    };
}


1;
