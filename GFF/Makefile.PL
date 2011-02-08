use strict;
use warnings;
use ExtUtils::MakeMaker;

WriteMakefile(
    NAME                => 'GFF',
    AUTHOR              => q{tnish <tnish@berkeley.edu>},
    VERSION_FROM        => 'lib/GFF.pm',
    ABSTRACT_FROM       => 'lib/GFF.pm',
    ($ExtUtils::MakeMaker::VERSION >= 6.3002
      ? ('LICENSE'=> 'perl')
      : ()),
    PL_FILES            => {},
    PREREQ_PM => {
        'Test::More' => 0,
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'GFF-*' },
);
