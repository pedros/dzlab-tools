#!/usr/bin/env perl
package main;
use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;
use File::Basename;
use PPI;
use List::Util 'max';
use Algorithm::Permute qw(permute);
use Statistics::Basic qw(:all nofill);

my $INH  = *ARGV;
my $ERRH = *STDERR;
my $OUTH = *STDOUT;
GetOptions(
    \%ARGV,
    'input|i=s',
    'output|o=s',
    'debug:i',
    'quiet'   => sub { $ARGV{quiet}   = 1; no diagnostics;  no warnings },
    'verbose' => sub { $ARGV{verbose} = 1; use diagnostics; use warnings },
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
) or pod2usage( -verbose => 1 );

IO:
{

    # We use the ARGV magical handle to read input
    # If user explicitly sets -i, put the argument in @ARGV
    if ( exists $ARGV{input} ) {
        unshift @ARGV, $ARGV{input};
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
}

my %data;
foreach my $file (@ARGV) {

    # Load a document from a file
    my $document = PPI::Document->new($file);

    # Strip out comments
    $document->prune('PPI::Token::Comment');

    # $document->normalized;
    # $document->index_locations;

    # Find all the named subroutines
    my $sub_nodes =
	$document->find( sub { $_[1]->isa('PPI::Statement::Sub') and $_[1]->name });

    next unless $sub_nodes;
    my %sub_names = map {
        $_->name => {
            'code'   => $_->content,
	    'tokens' => count($_)
	}
    } @$sub_nodes;

    my $file_name = fileparse $file;
    $data{$file_name} = \%sub_names;
}



#     foreach my $sub_name (keys %sub_names) {
# 	eval $sub_names{$sub_name};
# 	use B::Deparse;
# 	my $deparse = B::Deparse->new("-p", "-sC");
# 	my $body = $deparse->coderef2text(\&$sub_name);          #This only works like half the time. No idea why.




#Set up @file_names and @sub_names for easy reference
my @file_names = sort keys %data;
my @sub_names;
foreach my $sub_hash ( values %data ) {
    foreach my $sub_name ( keys %{$sub_hash} ) {
        push @sub_names, $sub_name;
    }
}



#remove duplicates and sort
my %hash = map { $_, 1 } @sub_names;
@sub_names = sort keys %hash;

#Put a second key into each subroutine hash, pointing to an array of the files the subroutine is used in.
#Each element in the array is a hash.
foreach my $sub_name (@sub_names) {
    my @file_names_by_sub;
    foreach my $file_name (@file_names) {
        if ( exists $data{$file_name}->{$sub_name} ) {
            push @file_names_by_sub, $file_name;
        }
    }
    foreach my $file_name_by_sub (@file_names_by_sub) {
        foreach my $file_name (@file_names) {
            $data{$file_name}->{$sub_name}->{'other files'}->{$file_name_by_sub} = $data{$file_name} 
              if exists $data{$file_name}->{$sub_name}
                  and $file_name ne $file_name_by_sub;   #keeps it from pointing to itself, and from creating new keys.
        }
    }
}

{
    local $Data::Dumper::Maxdepth = 4;
  #  print Dumper \%data;
 #   exit;
}



#cross correlation part 

my %cc_by_sub_name;

foreach my $sub_name (@sub_names) {
    my %code_of_sub_variations;
    foreach my $file ( keys %data ) {
        foreach my $my_sub_name ( keys %{ $data{$file} } ) {
            if ( $sub_name eq $my_sub_name) {
                $code_of_sub_variations{$file} =
                  $data{$file}->{$my_sub_name}->{'tokens'};
            }
        }
    }
    $cc_by_sub_name{$sub_name} = cc_all_combos(%code_of_sub_variations);
}



print Dumper \%cc_by_sub_name;



# Takes a hash representing a single subroutine. keys are files, values are hashes of the token frequencies.
# Returns a hash with keys being a filename pointing to a hash whose keys are another filename and whose value is the token comparison for the sub between those two files.
# So to retieve the difference between file2 and file3, do $hash{file2}{file3} or $hash{file3}{file2}.
sub cc_all_combos {
    my (%hash) = @_;
    my @files = keys %hash;
    my %ret_hash;
    my @combinations;

    while (my $file1 = shift @files) {
	foreach my $file2 (@files) {
	    my @combo = ($file1, $file2);
	    push @combinations, \@combo;
	}
    }

    for my $combo (@combinations) {
	my $file1 = @{$combo}[0];
	my $file2 = @{$combo}[1];
	$ret_hash{$file1}->{$file2} = cross_correlation($hash{$file1}, $hash{$file2});
	$ret_hash{$file2}->{$file1} = cross_correlation($hash{$file1}, $hash{$file2});
    }
    return \%ret_hash;
}



sub normalized_hd {
  my ($k, $l) = @_;
  my $xor = $k ^ $l;
  my $hamming_distance = $xor =~ tr/\0//c;
  my $max_length = max length $k, length $l;
  return sprintf ("%g", 1 - $hamming_distance / $max_length)
      unless $max_length == 0;
  return -1;
}


#P value? (list squares fit)
sub cross_correlation {
    my ($tokens1, $tokens2) = @_; 
    
    my %all_tokens = map { $_ => 0 } (keys %$tokens1, keys %$tokens2);
    my (@vec1, @vec2);

    for my $token (keys %all_tokens) { 

	$tokens1->{$token} = 0 unless $tokens1->{$token};
	$tokens2->{$token} = 0 unless $tokens2->{$token};

	push @vec1, $tokens1->{$token};
	push @vec2, $tokens2->{$token};
    }	

    my $cc = correlation(vector(@vec1), vector(@vec2));
    return "$cc";
}





#takes a sub_node object, returns a ref to a hash of token frequencies
sub count {
    my ($sub_node) = @_;
    
    my @tokens = $sub_node->find( 
	sub {
	    $_[1]->isa('PPI::Token')
		and not $_[1]->isa('PPI::Token::Whitespace')
	} );
    my %counts;
    $counts{ join "\t", ref ($_), $_->content }++ for @{$tokens[0]}; 
    return \%counts;    
}



sub norm {
    my ( $string, @strings ) = @_;
    my $longest = max( map { length } @strings );
    until ( length $string == $longest ) {
        $string = "$string ";
    }
    return $string;
}


sub print_table {
##print out a table showing which files contain which subroutines.
    my (@sub_names, @file_names, %data) = @_;
    print join( "\t", q{}, @sub_names ), "\n";
    
    foreach my $file_name (@file_names) {
	print $file_name, "\t";
	foreach my $sub_name (@sub_names) {
	    if ( $data{$file_name}->{$sub_name} ) {
		print q{*};
	    }
	    else { print q{}; }
	    print "\t";
	}
	print "\n";
    }
}



__DATA__

__END__
=head1 NAME

 APerlyName.pl - Short description

=head1 SYNOPSIS

 APerlyName.pl [OPTION]... [FILE]...

=head1 DESCRIPTION

 Long description

=head1 OPTIONS

 -i, --input       filename to read from                            (STDIN)
 -o, --output      filename to write to                             (STDOUT)
     --debug       print additional information
     --verbose     print diagnostic and warning messages
     --quiet       print no diagnostic or warning messages
     --version     print current version
     --license     print author's contact and copyright information
     --help        print this information
     --manual      print the plain old documentation page

=head1 VERSION

 0.0.1

=head1 REVISION

 $Rev: $:
 $Author: $:
 $Date: $:
 $HeadURL: $:
 $Id: $:

=head1 AUTHOR

 Pedro Silva <pedros@berkeley.edu/>
 Zilberman Lab <http://dzlab.pmb.berkeley.edu/>
 Plant and Microbial Biology Department
 College of Natural Resources
 University of California, Berkeley

=head1 COPYRIGHT

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut
