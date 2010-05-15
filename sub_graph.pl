#!/usr/bin/env perl

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

GetOptions(
    _meta_options (\%ARGV),
    # arguments go here
    # not above _meta_options!
) or pod2usage( -verbose => 1 );

my ( $ERRH, $INH, $OUTH ) = _prepare_io( \@ARGV );

my $data = load_files(\@ARGV);

my @file_names = extract_file_names ($data);
my @sub_names  = extract_sub_names  ($data);

#put pair_value key into sub hash.
foreach my $sub_name (@sub_names) {
    my @file_names_by_sub = _extract_files_by_sub ($sub_name, $data);
    _insert_pair_values($sub_name, @file_names_by_sub);
}
remove_token_keys($data);


{
    local $Data::Dumper::Maxdepth = 4;
#    die Dumper $data->{'extract_sequences.pl'}{'index_fasta'}{'pair value'}; 
    #print {$OUTH} Dumper $data;
    #exit;
}


auto_module( "/home/jgraff/workspace/bisulfite/trunk/module_test.pm" );  



  



##More stuff to do: add mode switches for table, subs of interest, and implement auto-generation of modules
#subs of interest thing should take a tolerance (0 - 1) and return all files/subs that exceed that number. Output: lists of subnames/filenames.
#auto-gen modules: find cases where correlation is perfect, and put those into a single module.
#Turn the whole thing into a module, make methods instead of modes.



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

sub load_files {
    my $ARGV = shift;
    
    my %data;
    foreach my $file (@$ARGV) {

        # Load a document from a file
        my $document = PPI::Document->new($file);

        # Strip out comments and documentation
        $document->prune('PPI::Token::Comment');
        $document->prune('PPI::Token::Pod');

        # Find all the named subroutines
        my $sub_nodes = $document->find(
            sub { $_[1]->isa('PPI::Statement::Sub') and $_[1]->name } );

        next unless $sub_nodes;

        my %sub_names = map {
            $_->name => {
                'code'   => $_->content,
                'tokens' => token_frequency($_)
            }
        } @$sub_nodes;

        my $file_name = fileparse $file;
        $data{$file_name} = \%sub_names;
    }

    return \%data;
}

sub extract_file_names {
    my ($data) = @_;

    return sort keys %$data;
}

sub extract_sub_names {
    my ($data) = @_;
    
    my @sub_names;
    foreach my $sub_hash ( values %$data ) {
        foreach my $sub_name ( keys %{$sub_hash} ) {
            push @sub_names, $sub_name;
        }
    }

    #remove duplicates and sort
    my %hash = map { $_, 1 } @sub_names;
    @sub_names = sort keys %hash;

    return @sub_names;
}


sub _extract_files_by_sub {
    my ($sub_name, $data) = @_;

    # all files that define current sub
    my @file_names_by_sub;
    foreach my $file_name (extract_file_names ($data) ) {
        if ( exists $data->{$file_name}->{$sub_name} ) {
            push @file_names_by_sub, $file_name;
        }
    }
    return @file_names_by_sub;
}

sub _insert_pair_values {
    my ($sub_name, @file_names_by_sub) = @_;
    foreach my $file_name_by_sub (@file_names_by_sub) {
	foreach my $file_name (extract_file_names ($data)) {
	    if (exists $data->{$file_name}->{$sub_name} and 
		$file_name ne $file_name_by_sub) {       
		my $cc = cross_correlation($data->{$file_name}->{$sub_name}->{'tokens'}, $data->{$file_name_by_sub}->{$sub_name}->{'tokens'});   #changed for cc
		$data->{$file_name}->{$sub_name}->{'pair value'}->{$file_name_by_sub}
		= $cc
	    }	    
	}
    }
}

sub remove_token_keys {
    my ($data) = @_;
    for my $file (@file_names) {
	for my $sub_name (keys %{$data->{$file}}) {
	    delete $data->{$file}{$sub_name}{'tokens'};
	}
    }
}


#I should examine the output more closely and make sure it matches up with the numbers.
sub auto_module {  
    my ( $module_file ) = @_;    
    my %cced_sub_code  = _get_perfect_cc_subs($data);
    unlink $module_file;
    _create_module( $module_file, sort keys %cced_sub_code )
        unless -e $module_file;
    _add_to_module( $cced_sub_code{$_}, $module_file ) foreach sort keys %cced_sub_code;
}
sub _get_perfect_cc_subs {     #This is returning multiple subs that are the same. problem? (I don't see how not to do that)
    my ($data) = @_;
    my %cced_sub_code;
    for my $sub_name (@sub_names) {    
	my $highest = 0;
	my @matched_files;
	for my $file1 (@file_names) {		
	    for my $file2 ( keys %{ $data->{$file1}{$sub_name}{'pair value'}}) {		
		my $cc = $data->{$file1}{$sub_name}{'pair value'}{$file2};
		unless ( $cc != 1 or scalar grep /^($file1|$file2)$/, @matched_files) {     #checks to make sure the cc is perfect, and neither of the two files has been matched previously for that sub.
		                                                                             #The grep thing helped a lot, but there are still a few doubles getting past.
		    push @matched_files, $file1, $file2;
		    my $code = $data->{$file1}{$sub_name}{'code'};		    
		    $code =~ s/./#$file1 \ns/; #comment specifying one file that it came from. If needed, the others can be found with this info.
		    if (exists $cced_sub_code{join '_', $sub_name, 1}) {			
			my $new_num = 1 + $highest;
			$code =~ s/sub ([^\s(]+)/sub $1_$new_num/;					
			$cced_sub_code{join '_', $sub_name, $new_num} = $code;
			$highest++;			
		    }		
		    elsif (exists $cced_sub_code{$sub_name}) { 		  #second time only
			my $new_code = delete $cced_sub_code{$sub_name}; 		    
			$new_code =~ s/sub ([^\s(]+)/sub $1_1/;		
			$cced_sub_code{join '_', $sub_name, 1} = $new_code;			                       
			redo;
		    } 
		    else {
			$cced_sub_code{$sub_name} = $code;
			$highest = 1;
		    }
		    last;
		}        
	    }		    		
	}   
    }
    return %cced_sub_code;
}

sub _add_to_module {
    my ( $code, $module ) = @_;
    open my $MODULE, '>>', $module or croak "Can't open $module: $!";
    print $MODULE "\n$code\n";
    close $MODULE or croak "Can't close $module: $!";
}
sub _create_module {
    my ( $file, @sub_names ) = @_;
    open my $FILE, '>', $file or croak "Can't open $file: $!";
    print $FILE
        "package Utils; \n",
        "use Export;\n",
        "our \@EXPORT_OK = qw(@sub_names);\n";    #basic module stuff

    close $FILE or croak "Can't close $file: $!";
}



#####This needs improvement! The groups are repeated w/ diff permutations. Possibly rethink approach. also broken now that subs_of_interest is gone.

#returns a hash, keys are sub names, values is an array of files.
#I can easily change this back to an array.
sub subs_of_interest {
    my ( $tolerance, $cc_by_sub_name, @sub_names ) = @_;
    my %cc_by_sub_name = %$cc_by_sub_name;
    my %subs_of_interest;
    for my $sub_name ( keys %cc_by_sub_name ) {
        my @groups;
        my @done;
        for my $file1 ( keys %{ $cc_by_sub_name{$sub_name} } ) {
            my @files = ($file1);
            for my $file2 ( keys %{ $cc_by_sub_name{$sub_name}{$file1} } ) {
                my $cc = $cc_by_sub_name{$sub_name}{$file1}{$file2};
                push @files, $file2
                    unless $cc < $tolerance
                        or grep /$file2/, @done;
                push @done, $file1;
            }
            push @groups, \@files
                if scalar @files > 1
            ; #if there's only one element, there are no matches for the given tolerance.
        }
        $subs_of_interest{$sub_name} = \@groups;

    }
    return %subs_of_interest;
}




sub cross_correlation {
    my ( $tokens1, $tokens2 ) = @_;

    my %all_tokens = map { $_ => 0 } ( keys %$tokens1, keys %$tokens2 );
    my ( @vec1, @vec2 );

    for my $token ( keys %all_tokens ) {

        $tokens1->{$token} = 0 unless $tokens1->{$token};
        $tokens2->{$token} = 0 unless $tokens2->{$token};

        push @vec1, $tokens1->{$token};
        push @vec2, $tokens2->{$token};
    }

    my $cc = correlation( vector(@vec1), vector(@vec2) );
    return "$cc";
}


#takes a sub_node object, returns a ref to a hash of token frequencies
sub token_frequency {
    my ($sub_node) = @_;

    my @tokens = $sub_node->find(
        sub {
            $_[1]->isa('PPI::Token')
                and not $_[1]->isa('PPI::Token::Whitespace');
        }
    );
    my %freqs;
    $freqs{ join "\t", ref($_), $_->content }++ for @{ $tokens[0] };
    return \%freqs;
}

sub print_table {
##print out a table showing which files contain which subroutines.
    my ( @sub_names, @file_names, $data ) = @_;
    print join( "\t", q{}, @sub_names ), "\n";

    foreach my $file_name (@file_names) {
        print $file_name, "\t";
        foreach my $sub_name (@sub_names) {
            if ( $data->{$file_name}->{$sub_name} ) {
                print q{*};
            }
            else { print q{}; }
            print "\t";
        }
        print "\n";
    }
}







#TRASH
#############################################################################

#OBSOLETE
sub norm {
    my ( $string, @strings ) = @_;
    my $longest = max( map {length} @strings );
    until ( length $string == $longest ) {
        $string = "$string ";
    }
    return $string;
}

#OBSOLETE
# Takes a hash representing a single subroutine. keys are files, values are hashes of the token frequencies.
# Returns a hash with keys being a filename pointing to a hash whose keys are another filename and whose value is the token comparison for the sub between those two files.
# So to retieve the difference between file2 and file3, do $hash{file2}{file3} or $hash{file3}{file2}.
sub cc_all_combos {
    my (%hash) = @_;
    my @files = keys %hash;
    my %ret_hash;
    my @combinations;

    while ( my $file1 = shift @files ) {
        foreach my $file2 (@files) {
            my @combo = ( $file1, $file2 );
            push @combinations, \@combo;
        }
    }

    for my $combo (@combinations) {
        my $file1 = @{$combo}[0];
        my $file2 = @{$combo}[1];
        $ret_hash{$file1}->{$file2}
            = cross_correlation( $hash{$file1}, $hash{$file2} );
        $ret_hash{$file2}->{$file1}
            = cross_correlation( $hash{$file1}, $hash{$file2} );
    }
    return \%ret_hash;
 }


#OBSOLETE
#Returns a hash representing a single subroutine. keys are files, values are hashes of the token frequencies.
sub get_sub_code_by_file {
    my ($sub_name, $data) = @_;
    my %code_of_sub_variations;
    foreach my $file ( keys %$data ) {
	foreach my $my_sub_name ( keys %{ $data->{$file} } ) {
	    if ( $sub_name eq $my_sub_name ) {
		$code_of_sub_variations{$file}
		= $data->{$file}->{$my_sub_name}->{'tokens'};
	    }
	}
    }
     return %code_of_sub_variations;
}


			#one way - add subs in sequentially
sub method1 {
    my ($cced_sub_code, $sub_name, $code, $highest) = @_;
    my $new_num = 1 + $highest;
    $code =~ s/sub ([^\s(]+)/sub $1_$new_num/;					
    $cced_sub_code->{join '_', $sub_name, $new_num} = $code;
}		       


			#another way - add a new sub as 1, and shift all the others up.
sub method2 {
    my ($cced_sub_code, $sub_name, $code) = @_;
    for my $sub_name_2 (reverse sort keys %$cced_sub_code) {
	if ($sub_name_2 =~ m/$sub_name/) {
	    my $new_sub_name = $sub_name_2;;
	    $new_sub_name =~ s/(\d+)$/$1+1/e;
	    $cced_sub_code->{$new_sub_name} = delete $cced_sub_code->{$sub_name_2};
	    $cced_sub_code->{$new_sub_name} =~ s/(\d+)([\s(])/join '', $1+1, $2/e;											     
	}
    }			
    $code =~ s/sub ([^\s(]+)/sub $1_1/;		
    $cced_sub_code->{$sub_name} = $code;	
}



#OBSOLETE
#cc part
# my %cc_by_sub_name;
# foreach my $sub_name (@sub_names) {
#     %code_of_sub_variations = get_sub_code_by_file($sub_name, $data);
#     my $sub_cc = cc_all_combos(%code_of_sub_variations);
#     insert_cc(%code_by_sub_variations, f
#     $data{$file1}{$sub_name}{'other files'}{$file2}{$cc}
#     #$cc_by_sub_name{$sub_name} = $sub_cc
#         if %$sub_cc;    #gets rid of subs used in only one file.
# }




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

 Jonathan Graff <jonthegraff@berkeley.edu>
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
