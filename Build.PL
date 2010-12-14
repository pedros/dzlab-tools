use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'DZLab::Tools',
    license             => 'gpl',
    dist_author         => q{Pedro Silva <pedros@berkeley.edu>},
    dist_version_from   => 'lib/DZLab/Tools.pm',
    build_requires => {
        'Test::More' => 0,
        'Test::Exception' => 0,
    },
    add_to_cleanup      => [ 'DZLab-Tools-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
