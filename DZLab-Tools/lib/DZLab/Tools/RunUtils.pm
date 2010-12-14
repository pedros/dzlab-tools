package DZLab::Tools::RunUtils;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Carp;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(_meta_options _prepare_io);

# use FindBin;
# use lib "$FindBin::Bin/DZLab-Tools/lib";
# use DZLab::Tools::RunUtils;

# * average_gff.pl
# * build_gff_score_table.pl
# * check_methyl.pl
# * eland2gff.pl
# * ensure_all_genes_present.pl
# * expand_gff.pl
# * gff2eland.pl
# * gff2eland3.pl
# * normalize_gff.pl
# * sort_gff.pl

=head1 EXPORTED FUNCTIONS

=head2 _meta_options $options

Takes a reference to a hash of boolean options with keys: quiet and verbose.

Returns a hash of commonly-used options for feeding into Getopt::Long.

=cut
sub _meta_options {
    my ($opt) = @_;

    return (
        'quiet'     => sub { $opt->{quiet}   = 1;          $opt->{verbose} = 0 },
        'verbose:i' => sub { $opt->{verbose} = $_[1] // 1; $opt->{quiet}   = 0 },
        'version'   => sub { pod2usage( -sections => ['VERSION', 'REVISION'],
                                        -verbose  => 99 )                      },
        'license'   => sub { pod2usage( -sections => ['AUTHOR', 'COPYRIGHT'],
                                        -verbose  => 99 )                      },
        'usage'     => sub { pod2usage( -sections => ['SYNOPSIS'],
                                        -verbose  => 99 )                      },
        'help'      => sub { pod2usage( -verbose  => 1  )                      },
        'manual'    => sub { pod2usage( -verbose  => 2  )                      },
    );
}

# * average_gff.pl
# * build_gff_score_table.pl
# * check_methyl.pl
# * eland2gff.pl
# * ensure_all_genes_present.pl
# * expand_gff.pl
# * gff2eland.pl
# * gff2eland3.pl
# * normalize_gff.pl
# parent_imprinting.pl
# * sort_gff.pl

=head2 _prepare_io $opt $argv

Takes a reference to a hash of options with keys: input, output, error, quiet; and a reference to an array of command line arguments.

Returns input, output and error file handles.

=cut
sub _prepare_io {
    my ($opt, $argv) = @_;

    my ($INH, $OUTH, $ERRH);
    
    unshift @$argv, $opt->{input} if exists $opt->{input};

    if (    exists $opt->{input} and exists $opt->{output}
               and $opt->{input} eq $opt->{output} ) {
        open $INH, q{<}, $opt->{input}
            or croak "Can't read $opt->{input}: $!";
        unlink $opt->{output};
    }
    else { $INH = *ARGV }

    if ( exists $opt->{output} and q{-} ne $opt->{output} ) {
        open $OUTH, q{>}, $opt->{output}
            or croak "Can't write $opt->{output}: $!";
    }
    else { $OUTH = *STDOUT }

    if ( exists $opt->{error} and q{-} ne $opt->{error} ) {
        open $ERRH, q{>}, $opt->{error}
            or croak "Can't write $opt->{error}: $!";
    }
    elsif ( exists $opt->{quiet} and $opt->{quiet} ) {
        use File::Spec;
        open $ERRH, q{>}, File::Spec->devnull
            or croak "Can't write $opt->{error}: $!";
    }
    else { $ERRH = *STDERR }

    return ( $INH, $OUTH, *STDERR = $ERRH );
}

1;


