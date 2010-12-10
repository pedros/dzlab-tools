#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use PPI;
use autodie;
use Digest::MD5 qw/md5_base64/;

my %subs = find_shared_subs(2,@ARGV);

foreach my $subname (keys %subs){
    my %md5;
    my @filenames = keys %{$subs{$subname}};

    # populate a hash of md5sums of content to [filenames]
    foreach my $filename (@filenames) {
        push @{$md5{$subs{$subname}{$filename}}}, $filename;
    }

    # for each md5sum (ie, identical subroutine definition),
    # print the sharing filenames and the actual duplicated subroutine
    # contents.
    say "#######################################################";
    say "#### $subname";
    foreach my $md5sum (keys %md5) {
        for my $f (sort @{$md5{$md5sum}}) {
            say "# $f";
        }
        
        say $subs{$subname}{$md5{$md5sum}->[0]};
        say "";
    }
}

# find_shared_subs($min_num_repeats, $file1, $file2, ...)
# 
# return hash of subnames to (hash of filenames to content)
# { subname => {file1 => content1, file2 => content2, ... }
sub find_shared_subs{
    my $min = shift;
    my %subs = ();

    for my $filename (@_){
        my $doc = PPI::Document->new($filename);
        $doc->prune('PPI::Token::Pod');
        $doc->prune('PPI::Token::Comment');

        $doc->find( 
            sub { 
                if ($_[1]->isa('PPI::Statement::Sub') and $_[1]->name){
                    my $subname = $_[1]->name;
                    my $content = $_[1]->content;

                    $subs{$subname}{$filename} = $content;
                }
            }
        );
    }
    return map {$_ => $subs{$_}} grep { keys %{$subs{$_}} >= $min } (keys %subs);
}

__END__

=head1 findsub.pl

 findsub.pl file1.pl file2.pl ...

=head1 Description

 Find copy/pasted subroutine definitions amongst perl files.

=cut
