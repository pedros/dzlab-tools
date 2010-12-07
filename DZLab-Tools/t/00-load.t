#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'DZLab::Tools' ) || print "Bail out!
";
}

diag( "Testing DZLab::Tools $DZLab::Tools::VERSION, Perl $], $^X" );
