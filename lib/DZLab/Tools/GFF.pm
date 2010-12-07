package DZLab::Tools::GFF;

use strict; use warnings;



sub gff_read {
    return [] if $_[0] =~ m/^
                            \s*
                            \#+
                           /mx;

    my ($seqname, $source, $feature, $start, $end,
        $score,   $strand, $frame,   $attribute
    ) = split m/\t/xm, shift || return;

    $attribute =~ s/[\r\n]//mxg;

    return {
        seqname   => lc $seqname,
        source    => $source,
        feature   => $feature,
        start     => $start,
        end       => $end,
        score     => $score,
        strand    => $strand,
        frame     => $frame,
        attribute => $attribute
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
