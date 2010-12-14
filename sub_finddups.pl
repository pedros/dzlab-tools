#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use PPI;
use autodie;

my %subs = find_shared_subs(2,@ARGV);

foreach my $subname (keys %subs){

    my @filenames = keys %{$subs{$subname}};

    # {file1 => content1, file2 => content2, ... }
    my %content = %{$subs{$subname}};

    # {file1 => normalized1, file2 => normalized2, ... }
    my %normalized = map { 
        $_ => PPI::Document->new(\$content{$_})->normalized 
    } @filenames;


    # group files into identical normalized forms.  We have to do this because
    # apparently normalized form only allows direct comparison with '=='...?
    my @groups;
    F:
    for my $i (0 .. $#filenames) {
        my $fi = $filenames[$i];

        # if the file is already found in one of the groups
        next F if (0 < grep { $_ eq $fi} 
            (map {@$_} @groups)); # flatten list

        my @g = ($fi);
        for my $j ($i + 1  .. $#filenames) {
            my $fj = $filenames[$j];

            # if normalized form equal, consider the two sub defs equal.
            if ($normalized{$fi} == $normalized{$fj}){
                push @g,$fj;
            }
        }
        push @groups, \@g;
    }

    say "##########################################################";
    say "### $subname";
    say "";

    for my $g (@groups) {
        for my $file (@$g) {
            say "# $file";
        }
        say $content{$g->[0]};
    }
    say "";
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
