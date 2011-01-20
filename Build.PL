use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'DZLab',
    license             => 'gpl3',
    dist_author         => [q{Pedro Silva <pedros@berkeley.edu>},
                            q{Toshiro Nishimura <tnish@berkeley.edu>}],
    dist_version        => '0.0.1',
    dist_abstract       => 'UC Berkeley Zilberman Lab Tools',
    requires => {
        'Algorithm::Permute'           => 0,
        # don't want to deal with these yet...
        #'Bio::DB::Fasta'               => 0,
        #'Bio::DB::Sam'                 => 0,
        #'Bio::SeqIO'                   => 0,
        'DBI'                          => 0,
        'DBD::SQLite'                  => 0,
        'File::Slurp'                  => 0,
        'Getopt::Long::Descriptive'    => 0,
        'List::MoreUtils'              => 0,
        'LWP::Simple'                  => 0,
        'Math::Round'                  => 0,
        'Moose'                        => 0,
        'Parallel::ForkManager'        => 0,
        'Pod::Usage'                   => 0,
        'Regexp::Common'               => 0,
        'Smart::Comments'              => 0,
        'Statistics::Basic'            => 0,
        'Statistics::KernelEstimation' => 0,
        #'System::Wrapper'              => 0,
        #'System::Wrapper::Parallel'    => 0,
        'Test::Exception'              => 0,
        'Test::More'                   => 0,
        #only for sub_graph.pm/pl, not in use
        #'XML::Simple'                  => 0,
        #'XML::Twig'                    => 0,
        'version'                      => 0,
    },
    script_files        => [
        glob("*.pl *.sh"), 
        "DZLab-Tools/gff2sqlite.pl"
    ],
    pm_files            => {
        'SequenceLength.pm' => 'lib/SequenceLength.pm',
        map { 
            my $original = $_; 
            s{^DZLab-Tools/}{}; # strip off DZLab-Tools/
            $original => $_
        } glob("DZLab-Tools/lib/DZLab/Tools/*.pm"),
    },
    add_to_cleanup      => [ 'DZLab-*' ],
);

$builder->create_build_script();

