use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'System::Wrapper',
    license             => 'gpl3',
    dist_author         => q{Pedro Silva <pedros@berkeley.edu>},
    dist_version_from   => 'lib/System/Wrapper.pm',
    build_requires => {
        'Test::More' => 0,
    },
    add_to_cleanup      => [ 'System-Wrapper-*' ],
    create_makefile_pl => 'traditional',
);

$builder->create_build_script();
