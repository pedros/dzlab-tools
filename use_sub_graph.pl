#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;
use lib q{/home/jgraff/workspace/bisulfite/trunk/Devel/Sub}; 
use Devel::Sub::Refactored;

#grab input from command line 


GetOptions(
    _meta_options (\%ARGV),
    # arguments go here
    # not above _meta_options!
) or pod2usage( -verbose => 1 );

my ( $ERRH, $INH, $OUTH ) = _prepare_io( \@ARGV );

my @files = @ARGV;



my $obj = Devel::Sub::Refactored->new('files' => \@files);
print $obj->load->pair_value( 'appender', 'cat_and_window.pl', 'test.pl');                               #auto_module("/home/jgraff/workspace/bisulfite/trunk/module_test.pm");  

#pair_value( 'appender', 'cat_and_window.pl', 'test.pl');


#auto_module("/home/jgraff/workspace/bisulfite/trunk/module_test.pm");  

# * = tested (not thoroughly though)


# print_table   *
# auto_module   *
# subs_of_interest/perfect_cc_subs  *
# show_cc_by_sub  * 


# selectors: 
# code($file_name, $sub_name)
# tokens($file_name, $sub_name)
# pair_value($sub_name, $file1, $file2)
# subs($file_name)
# files($sub_name)

    
    
    
#selectors/mutators for the stuff in $fields






sub _meta_options {
    my $OPT = @_;
    
    return (
        $OPT,
        'input|i=s',
        'output|o=s',
        'debug:i',
        'quiet'   => sub { $OPT->{quiet}   = 1; no diagnostics;  no warnings },
        'verbose' => sub { $OPT->{verbose} = 1; use diagnostics; use warnings },
        'version' => sub {
            pod2usage
                -sections => [ 'VERSION', 'REVISION' ],
                    -verbose  => 99;
        },
        'license' => sub {
            pod2usage
                -sections => [ 'AUTHOR', 'COPYRIGHT' ],
                    -verbose  => 99;
        },
        'usage' => sub {
            pod2usage
                -sections => ['SYNOPSIS'],
                    -verbose  => 99;
        },
        'help'   => sub { pod2usage -verbose => 1 },
        'manual' => sub { pod2usage -verbose => 2 },
    );
}

sub _prepare_io {
    my $ARGV = shift;

    my $INH  = *ARGV;
    my $ERRH = *STDERR;
    my $OUTH = *STDOUT;

    # We use the ARGV magical handle to read input
    # If user explicitly sets -i, put the argument in @$ARGV
    if ( exists $ARGV{input} ) {
        unshift @$ARGV, $ARGV{input};
    }

    # Allow in-situ arguments (equal input and output filenames)
    # FIXME: infinite loop. Why?
    if (    exists $ARGV{input}
        and exists $ARGV{output}
        and $ARGV{input} eq $ARGV{output} )
    {
        croak "Bug: don't use in-situ editing (same input and output files";
        open $INH, q{<}, $ARGV{input}
            or croak "Can't read $ARGV{input}: $!";
        unlink $ARGV{input};
    }

    # Redirect STDOUT to a file if so specified
    if ( exists $ARGV{output} and q{-} ne $ARGV{output} ) {
        open $OUTH, q{>}, $ARGV{output}
            or croak "Can't write $ARGV{output}: $!";
    }

    return ( $ERRH, $INH, $OUTH );
}
