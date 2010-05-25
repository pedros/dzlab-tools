#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use version; our $VERSION = qv('0.0.1');

use File::Path qw(make_path remove_tree);
use File::Spec;

use threads;

GetOptions(
    \%ARGV,
    'input|i=s',        'output|o=s',       'error|e=s',
    'ecotype-a|ea=s',   'ecotype-b|eb=s',   'genotype|g=s', 'tissue|t=s',
    'reference-a|ra=s', 'reference-b|rb=s', 'annotation|a=s',
    _meta_options( \%ARGV ),
)
and (@ARGV or $ARGV{input})
and @ARGV{ qw/
                 ecotype-a ecotype-b genotype tissue
                 reference-a reference-b annotation
             / }
or pod2usage( -verbose => 1 );

my ( $INH, $OUTH, $ERRH ) = _prepare_io( \%ARGV, \@ARGV );

my %file_names = (
    a => build_file_names(
        @ARGV{ qw/input output ecotype-a ecotype-b tissue genotype/ } ),
    b => build_file_names(
        @ARGV{ qw/input output ecotype-b ecotype-a tissue genotype/ } ),
);

my %commands = (
    a => build_common_commands( 
        $file_names{a}, @ARGV{ qw/reference-a annotation/ } ),
    b => build_common_commands( 
        $file_names{b}, @ARGV{ qw/reference-b annotation/ } ),
    c => build_unique_commands(
        $file_names{a}, $file_names{b}, @ARGV{ qw/ecotype-a ecotype-b/ }),
);

@{$commands{a}} == @{$commands{b}} and @{$commands{a}} == @{$commands{c}}
or croak "Wrong number of job stages";

build_output_directory( @ARGV{ qw/output ecotype-a ecotype-b/ } );

while ( @{$commands{a}} and @{$commands{b}} and @{$commands{c}} ) {

    my %threads = (
        a => threads->new( \&run_commands, shift @{$commands{a}} ),
        b => threads->new( \&run_commands, shift @{$commands{b}} ),
    );

    $threads{a}->join;
    $threads{b}->join;

    $threads{c} = threads->new( \&run_commands, shift @{$commands{c}} );
    $threads{c}->join;
}


sub build_output_directory {
    my ($out_dir) = @_;

    if (-d $out_dir
        and print STDERR "Overwrite $out_dir (y/N)?: " and q{y} eq <STDIN>) {
        ### removing $out_dir
        remove_tree ( $out_dir, {keep_root => 1} );
    }
    else {
        make_path ( File::Spec->catfile ($out_dir) );
    }
}

sub build_unique_commands {
    my ($names_a, $names_b) = @_;

    my @pre_wrapup = (

        "split_on_mismatches.pl -a $ARGV{'ecotype-a'} -b $ARGV{'ecotype-b'} \\
             $names_a->{eland} $names_b->{eland}",
    );

    my @post_wrapup = (

    );
    
    return [\@pre_wrapup, \@post_wrapup];
}

sub build_common_commands {
    my ($names, $reference, $annotation) = @_;

    my @pre_processing = (

        need_file( $reference . '.1.ebwt' )
        ? "bowtie-build --quiet $reference $reference" : 0,
    
        need_file( $names->{bowtie} )
        ? "bowtie -B 1 --quiet -v 3 --best \\
            $reference $names->{fastq} $names->{bowtie}" : 0,

        need_file( $names->{eland} )
        ? "parse_bowtie.pl -u $names->{fastq} -o $names->{eland}  \\
            $names->{bowtie}       " : 0,

        "sort -k1,1 -S 50% $names->{eland} > $names->{eland}.sort \\
            && mv $names->{eland}.sort $names->{eland}",
    );


    my @post_processing = (

        need_file( $names->{reads_gff} )
        ? "parse_eland.pl -3 -o $names->{reads_gff} $names->{eland} " : 0,

        "sort -k1,1 -k4,4n -k5,5n -k7,7 -S 50% \\
            $names->{reads_gff} > $names->{reads_gff}.sort \\
            && mv $names->{reads_gff}.sort $names->{reads_gff}",

        need_file( $names->{repeats} )
        ? "filter_repeats.pl -o $names->{reads_gff}.filter \\
            $names->{reads_gff} 2> $names->{repeats} \\
            && mv $names->{reads_gff}.filter $names->{reads_gff}" : 0,

        need_file( $names->{genes_gff} )
        ? "window_gff_REFACTORED.pl -c sum -k -t ID -f gene -r \\
            -g $annotation -o $names->{genes_gff} $names->{reads_gff}" : 0,
    );

    return [\@pre_processing, \@post_processing];
         
}

sub build_file_names {
    my ($input, $out_dir, $ecotype_a, $ecotype_b, $tissue, $genotype) = @_;

    my $base_name = File::Spec->catfile(
        $out_dir,
        join (q{_},
              $ecotype_a . q{x} . $ecotype_b,
              $tissue,
              $genotype,
          )
    );

    my $fastq_name   = $input;
    my $bowtie_names = join q{_}, $base_name, $ecotype_a . '.bowtie';
    my $eland_names  = join q{_}, $base_name, $ecotype_a . '.eland3';
    my $gff_names    = join q{_}, $base_name, $ecotype_a . '.gff';
    my $repeat_names = join q{_}, $base_name, $ecotype_a . '.repeats';
    my $gff_genes    = join q{_}, $base_name, $ecotype_a, 'genes' . '.gff';

    return {
        fastq     => $fastq_name,
        bowtie    => $bowtie_names,
        eland     => $eland_names,
        reads_gff => $gff_names,
        genes_gff => $gff_genes,
        repeats   => $repeat_names,
    };
}

sub need_file {
    !-f $_[0] or !-s $_[0];
}

sub run_commands {
    my ($commands) = @_;

    run_cmd( $_ ) for grep { $_ } @$commands;

    return 1;
}


sub run_cmd {
    my ($cmd) = @_;

    print STDERR "-- CMD: $cmd\n";

    my $exit_code = system $cmd;
    croak "** FAIL: non-zero exit status ($exit_code)" if $exit_code;

    return 1;
}


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

sub _prepare_io {
    my ($opt, $argv) = @_;

    my ($INH, $OUTH, $ERRH);
    
    # If user explicitly sets -i, put the argument in @$argv
    unshift @$argv, $opt->{input} if exists $opt->{input};

    # Allow in-situ arguments (equal input and output filenames)
    if (    exists $opt->{input} and exists $opt->{output}
               and $opt->{input} eq $opt->{output} ) {
        open $INH, q{<}, $opt->{input}
            or croak "Can't read $opt->{input}: $!";
        unlink $opt->{output};
    }
    else { $INH = *ARGV }

    $OUTH = *STDOUT;
    my @date = localtime (time);
    $opt->{output} //= File::Spec->catdir (
        File::Spec->curdir,
        sprintf "DZ_full-run_%4d-%02d-%02d_%02d.%02d.%02d",
        $date[5] + 1900, $date[4] + 1, $date[3], $date[2], $date[1], $date[0]
    );

    # Log STDERR if so specified
    if ( exists $opt->{error} and q{-} ne $opt->{error} ) {
        open $ERRH, q{>}, $opt->{error}
            or croak "Can't write $opt->{error}: $!";
    }
    elsif ( exists $opt->{quiet} and $opt->{quiet} ) {
        open $ERRH, q{>}, File::Spec->devnull
            or croak "Can't write $opt->{error}: $!";
    }
    else { $ERRH = *STDERR }

    return ( $INH, $OUTH, *STDERR = $ERRH );
}

__DATA__


__END__

=head1 NAME

 APerlyName.pl - Short description

=head1 SYNOPSIS

 APerlyName.pl [OPTION]... [[-i] FILE]...

=head1 DESCRIPTION

 Long description

=head1 OPTIONS

 -i, --input       <string>     input filename                           (STDIN)
 -o, --output      <string>     output filename                          (STDOUT)
 -e, --error       <string>     output error filename                    (STDERR)
     --verbose     [integer]    print increasingly verbose error messages
     --quiet                    print no diagnostic or warning messages
     --version                  print current version
     --license                  print author's contact and copyright information
     --help                     print this information
     --manual                   print the plain old documentation page

=head1 VERSION

 0.0.1

=head1 REVISION

 $Rev: $:
 $Author: $:
 $Date: $:
 $HeadURL: $:
 $Id: $:

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
