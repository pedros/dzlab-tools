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
use Digest::MD5 qw(md5_hex);
use System::Wrapper;
use System::Wrapper::Parallel;

use threads;

GetOptions(
    \%ARGV,
    'input|i=s',        'output|o=s',       'error|e=s',
    'ecotype-a|ea=s',   'ecotype-b|eb=s',   'genotype-a|ga=s', 'genotype-b|gb=s', 
    'tissue|t=s', 'loci-type|l=s',
    'reference-a|ra=s', 'reference-b|rb=s', 'annotation|a=s', 'debug',
    'overwrite',
    _meta_options( \%ARGV ),
)
and (@ARGV or $ARGV{input})
and @ARGV{ qw/
                 ecotype-a ecotype-b genotype tissue
                 reference-a reference-b annotation
             / }
or pod2usage( -verbose => 1 );

$ARGV{overwrite} //= 0;

my ( $INH, $OUTH, $ERRH ) = _prepare_io( \%ARGV, \@ARGV );

my %file_names = (
    a => build_file_names(
        @ARGV{ qw/input output ecotype-a ecotype-b tissue genotype-a genotype-b loci-type/ }, 0 ),
    b => build_file_names(
        @ARGV{ qw/input output ecotype-a ecotype-b tissue genotype-a genotype-b loci-type/ }, 1 ),
);

my %commands = (
    a => build_common_commands( $file_names{a}, @ARGV{ qw/reference-a annotation loci-type/ } ),
    b => build_common_commands( $file_names{b}, @ARGV{ qw/reference-b annotation loci-type/ } ),
    c => build_unique_commands( $file_names{a}, $file_names{b} ),
);

my $done_dir 
= build_output_directory( @ARGV{ qw/output ecotype-a ecotype-b/ } );

while ( @{$commands{a}} ) {

    my %threads = (
        a => threads->new( \&run_commands, shift @{$commands{a}}, $done_dir ),
        b => threads->new( \&run_commands, shift @{$commands{b}}, $done_dir ),
    );

    $threads{a}->join;
    $threads{b}->join;

    $threads{c} = threads->new( \&run_commands, shift @{$commands{c}}, $done_dir );
    $threads{c}->join;
}



sub build_output_directory {
    my ($out_dir) = @_;

    my $done_dir
    = File::Spec->catfile (File::Spec->catdir($out_dir, '.done'));

    if (-d $out_dir and $ARGV{overwrite}) {
        remove_tree ( $out_dir, {keep_root => 1} );
    }
    else {
        make_path ( $done_dir );
    }
    return $done_dir;
}

sub build_unique_commands {
    my ($names_a, $names_b) = @_;

    my @pre_wrapup = (

        (
            need_file( $names_a->{reads_gff} )
         || need_file( $names_b->{reads_gff} )
        )
        ? System::Wrapper->new(
            interpreter => 'perl',
            executable  => 'split_on_mismatches.pl',
            arguments   => [ -a => $ARGV{'ecotype-a'}, -b => $ARGV{'ecotype-b'},
                             $names_a->{eland}, $names_b->{eland} ])
        : 0,
    );
    
    my @post_wrapup = (
        need_file( $names_a->{table} )
        ? System::Wrapper->new(
            interpreter => 'perl',
            executable  => 'build_gff_score_table.pl',
            input       => [
                $names_a->{genes_gff},
                $names_a->{genes_filter},
                $names_b->{genes_gff},
                $names_b->{genes_filter},
            ],
            output      => { -o => $names_a->{table} },)
        : 0,
    );
    
    return [\@pre_wrapup, \@post_wrapup];
}

sub build_common_commands {
    my ($names, $reference, $annotation, $loci_type) = @_;

    $loci_type ||= 'genes';
    $loci_type =~ s/s$//;

    my @pre_processing = (

        need_file( $reference . '.1.ebwt' )
        ? System::Wrapper->new(
            executable => 'bowtie-build',
            arguments  => ['--quiet', $reference, $reference],
        )
        : 0,
    
        need_file( $names->{bowtie} )
        ? System::Wrapper->new(
            executable => 'bowtie',
            arguments => [ -B => 1, '--quiet', -v => 3, '--best' ],
            input     => [$reference, $names->{fastq}],
            output    => [$names->{bowtie}])
        : 0,

        need_file( $names->{eland} )
        ? System::Wrapper->new(
            interpreter => 'perl',
            executable  => 'parse_bowtie.pl',
            arguments   => [ -u => $names->{fastq} ],
            input       => [ $names->{bowtie} ],
            output      => { -o => $names->{eland}})
        : 0,

        need_file( $names->{reads_gff} )
        ? System::Wrapper->new(
            executable => 'sort',
            arguments  => [ -k => '1,1', -S => '15%'],
            input      => [ $names->{eland} ],
            output     => { q{>} => $names->{eland} },)
        : 0,
    );

    my @post_processing = (

        need_file( $names->{reads_gff} )
        ? System::Wrapper->new(
            interpreter => 'perl', 
            executable  => 'parse_eland.pl',
            arguments   => [ '-3' ],
            input       => [ $names->{eland} ],
            output      => { -o => $names->{reads_gff}},)
        : 0,

        need_file( $names->{repeats} )
        ? System::Wrapper->new(
            executable => 'sort',
            arguments  => [ -k => '1,1', -k => '4,4n', -k => '5,5n', -k => '7,7', -S => '15%' ],
            input      => [ $names->{reads_gff}],
            output     => { q{>} => $names->{reads_gff} },)
        : 0,

        need_file( $names->{filter_gff} )
        ? System::Wrapper->new(
            interpreter => 'perl',
            executable  => 'filter_repeats.pl',
            output      => {
                -o   => $names->{filter_gff},
                '2>' => $names->{repeats},
            },
            input       => [ $names->{reads_gff} ],)
        : 0,

        need_file( $names->{genes_gff} )
        ? System::Wrapper->new(
            interpreter => 'perl',
            executable  => 'window_gff.pl',
            arguments   => [ -c => 'sum', '-k', -t => 'ID', -f => $loci_type, '-r',
                             -g => $annotation ],
            output      => { -o => $names->{genes_gff} },
            input       => [ $names->{reads_gff} ],)
        : 0,

        need_file( $names->{genes_filter} )
        ? System::Wrapper->new(
            interpreter => 'perl',
            executable  => 'window_gff.pl',
            arguments   => [ -c => 'sum', '-k', -t => 'ID', -f => $loci_type, '-r',
                             -g => $annotation ],
            output      => { -o => $names->{genes_filter} },
            input       => [ $names->{filter_gff} ],)
        : 0,
    );

    return [\@pre_processing, \@post_processing];
         
}

sub build_file_names {
    my ($input, $out_dir, $ecotype_a, $ecotype_b, $tissue, $genotype_a, $genotype_b, $loci_type, $reverse) = @_;

    my $base_name = File::Spec->catfile(
        $out_dir,
        join (q{_},
              "$ecotype_a-$genotype_a" . q{x} . "$ecotype_b-$genotype_b",
              $tissue,
          )
    );

    $loci_type ||= 'genes';
    my $ecotype = $reverse ? $ecotype_b : $ecotype_a;

    my $fastq_name   = $input;
    my $bowtie_names = join q{_}, $base_name, $ecotype . '.bowtie';
    my $eland_names  = join q{_}, $base_name, $ecotype . '.eland3';
    my $gff_names    = join q{_}, $base_name, $ecotype . '.gff';
    my $filter_names = join q{_}, $base_name, $ecotype . '_filtered.gff';
    my $gff_genes    = join q{_}, $base_name, $ecotype . ".$loci_type.gff";
    my $filter_genes = join q{_}, $base_name, $ecotype . "_filtered.$loci_type.gff";
    my $repeat_names = join q{_}, $base_name, $ecotype . ".$loci_type.repeats";
    my $score_table  = join q{_}, $base_name           . ".$loci_type.table";

    return {
        fastq        => $fastq_name,
        bowtie       => $bowtie_names,
        eland        => $eland_names,
        reads_gff    => $gff_names,
        filter_gff   => $filter_names,
        genes_gff    => $gff_genes,
        genes_filter => $filter_genes,
        repeats      => $repeat_names,
        table        => $score_table,
    };
}

sub need_file {
    !-f $_[0] or !-s $_[0];
}

sub run_commands {
    my ($commands, $done_dir) = @_;

    return unless ref $commands and 'ARRAY' eq ref $commands;
    
    for (grep { $_ } @$commands) {
        $_->debug   = 1 if $ARGV{debug};
        $_->verbose = 1 if $ARGV{verbose};

        my $out_md5 = File::Spec->catfile( $done_dir, md5_hex( "$_ " ) );

        unless (-e $out_md5) {
            my $stdout = $_->run;

            if (defined ($stdout) and 0 == $stdout ) {
                open my $DONE, q{>}, $out_md5 or die $!;
                print $DONE "$_";
            }
        }
    }

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

 -i   --input         short reads file in fastq format
 -o   --output        output filename
 -e   --error         output error filename 
 -ea  --ecotype-a     label for maternal strain (eg. columnbia)
 -eb  --ecotype-b     label for paternal strain (eg. ler)
 -ga  --genotype-a    label for ecotype-a genotype (WT, mut-a, mut-b, ...)
 -gb  --genotype-b    label for ecotype-b genotype (WT, mut-a, mut-b, ...)
 -t   --tissue        label for tissue type (endosperm, embryo, etc)
 -ra  --reference-a   genome of ecotype-a in fasta format
 -rb  --reference-b   genome of ecotype-b in fasta
 -a   --annotation    gff annotation 
 -l   --loci-type     [gene, transposon, exon, etc]

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
