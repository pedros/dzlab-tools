use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'DZLab::Tools',
    license             => 'gpl3',
    dist_author         => [q{Pedro Silva <pedros@berkeley.edu>},
                            q{Toshiro Nishimura <tnish@berkeley.edu>}],
    dist_version_from   => 'lib/DZLab/Tools.pm',
    build_requires => {
        'Test::More' => 0,
        'Test::Exception' => 0,
    },
    requires       => {
        'strict'     => 0,
        'warnings'   => 0,
        'version'    => 0,
        'Carp'       => 0,
    },
    add_to_cleanup      => [ 'DZLab-Tools-*' ],
    create_makefile_pl  => 'traditional',
    sign                => 1,
    create_readme       => 1,
    create_license      => 1,
);

$builder->create_build_script();
