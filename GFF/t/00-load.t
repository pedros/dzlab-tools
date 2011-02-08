#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'GFF' ) || print "Bail out!
";
}

diag( "Testing GFF $GFF::VERSION, Perl $], $^X" );
