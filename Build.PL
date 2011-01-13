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
        'Test::More' => 0,
        'version'    => 0,
    },
    script_files        => [glob("*.pl *.sh")],
    add_to_cleanup      => [ 'DZLab-*' ],
);

$builder->create_build_script();
