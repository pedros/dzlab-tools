#!/usr/bin/env perl
package Devel::Sub::Refactored;

use warnings;
use strict;
use Getopt::Long qw(:config gnu_getopt);
use Pod::Usage;
use File::Basename;
use PPI;
use Statistics::Basic qw(:all nofill);
use Carp;
use Data::Dumper;

#verbose option stuff?

my $fields = {
    keep_comments => 0,
    files         => [],
    subs          => [],
    pair_values   => 1,
    keep_tokens   => 0,
};


#### data-handling methods ####


#Only initializes the object. load everything with the load method.
#files can be passed here, or later to load.
sub new {
    my ($class, %args) = @_;
    my $self = bless $fields, ref $class ? ref $class : $class;
    while (my ($key, $value) = each %args) {
	croak "Can't access '$key' field in class $class"
	    if !exists $fields->{$key}
	or $key =~ m/^_/;
	$self->{$key} = $value;
    }
    return $self;
}


# uses obj->files unless it is passed files as args.
sub load {
    my ($self, @files) = @_;    
    return 
	unless @files or ($self->files and 'ARRAY' eq ref $self->files);
    
    $self->_set_files_instance_var(\@files);

    @files = $self->files;

    $self->_load_files(@files);
    
    my $data = $self->{data};
    $self->_set_subs_instance_var(_extract_sub_names($data) );   # for easy reference
    
    $self->insert_pair_values
	unless $self->no_pair_values;
    
    $self->_remove_token_keys
	unless $self->keep_tokens;

    return $self;
}


#should I call this automatically?
#put pair_value key into sub hash.
sub insert_pair_values {
    my ($self) = @_;
    foreach my $sub_name ($self->subs) {
	my @file_names_by_sub = $self->files($sub_name);
	foreach my $file_name_by_sub (@file_names_by_sub) {
	    foreach my $file_name ($self->files) {
		if (exists $self->{data}{$file_name}->{$sub_name} and 
		    $file_name ne $file_name_by_sub) {       
		    my $cc = _cross_correlation($self->{data}{$file_name}{$sub_name}{'tokens'}, $self->{data}{$file_name_by_sub}{$sub_name}{'tokens'});  
		    $self->{data}{$file_name}{$sub_name}{'pair value'}{$file_name_by_sub}
		    = $cc
		}	    
	    }
	}
    }
    return $self;
}



#### OUTPUT METHODS ####

#I should examine the output more closely and make sure it matches up with the numbers.
sub auto_module {  
    my ( $self, $module_file ) = @_;    
    my %cced_sub_code  = _get_perfect_sub_code($self);
    unlink $module_file;  #for convienience
    _create_module( $module_file, sort keys %cced_sub_code )
        unless -e $module_file;
    _add_to_module( $cced_sub_code{$_}, $module_file ) foreach sort keys %cced_sub_code;
}


sub print_table {
##print out a table showing which files contain which subroutines.
    my ( $self ) = @_;
    print join( "\t", q{}, $self->subs ), "\n";

    foreach my $file_name ($self->files) {
        print $file_name, "\t";
        foreach my $sub_name ($self->subs) {
            if ( $self->{data}{$file_name}->{$sub_name} ) {
                print q{*};
            }
            else { print q{}; }
            print "\t";
        }
        print "\n";
    }
}


sub show_cc_by_sub {
    my ($self) = @_;
    my %cc_by_sub;

    for my $file1 ($self->files) {
	for my $sub_name (keys %{$self->{data}{$file1}}) {
	    for my $file2 (keys %{$self->{data}{$file1}{$sub_name}{'pair value'}} ) {    #should be (keys $self->pair_value($sub_name, $file1)), but that doesn't work. why?
		my $cc = $self->pair_value($sub_name, $file1, $file2);
		$cc_by_sub{$sub_name}{$file1}{$file2} = $cc;
	    }}}
    print Dumper \%cc_by_sub;
    exit;
}



sub get_perfect_cc_subs {       #Only works for a tolerance of 1. In order to make it work for others, it would have to return many permutations, which we don't want. 
                                 #Perhaps if I didn't do it by groups, but just did it one-to-one? But then it wouldn't be compatible with get_perfect_sub_code.
    my ($self) = @_;
    my %cced_sub_code;
    my %subs_of_interest;
    for my $sub_name ($self->subs) {    
	my $highest = 0;
	my @matched_files;
	my @groups;
	for my $file1 ($self->files) {	
	    my @files = ($file1);
	    for my $file2 (keys %{$self->{data}{$file1}{$sub_name}{'pair value'}} ) {   #should be (keys $self->pair_value($sub_name, $file1)), but that doesn't work. why?
		my $cc = $self->pair_value($sub_name, $file1, $file2);
		if ( $cc == 1 and not grep /^($file1|$file2)$/, @matched_files) {     #checks to make sure the cc is high enough, and neither of the two files has been matched previously for that sub.
		    push @files, $file2;
		    push @matched_files, $file2;		    
		}		
	    }
	    push @groups, \@files
		if scalar @files > 1;
	}
	$subs_of_interest{$sub_name} = \@groups
	    if scalar @groups > 0;
    }   
    return %subs_of_interest;
}



#### Selectors / Mutators ####

sub code {
    my ($self, $file, $sub_) = @_;
    return $self->{data}{$file}{$sub_}{'code'};
}


#if not given a second file, returns the whole pair value hash for the given file and sub.
sub pair_value {
    my ( $self, $sub_, $file1, $file2) = @_;
    return $self->{data}{$file1}{$sub_}{'pair value'}{$file2} if defined $file2;
    
    return $self->{data}{$file1}{$sub_}{'pair value'};
}

#if given files, rets all subs defined in those files. If given no files, rets all subs defined in all files.
sub subs {
    my ($self, @files) = @_;
    
    unless (@files) {

	return wantarray ? @{$self->{subs}} : $self->{subs};
    }
    my @subs;
    for my $file (@files) {
	push @subs, sort keys %{ $self->{data}{$file}};
    }
    return @subs;
}


#rets all files the given subs are defined in. If given no subs, rets all files.
sub files {
    my ($self, @subs) = @_;

    unless (@subs) {
	return wantarray ? @{$self->{files}} : $self->{files};
    }

    my @files;
    for my $sub_ (@subs) {	
	for my $file (keys %{$self->{data}}) {	    
	    push @files, $file
		if exists $self->{data}{$file}{$sub_};        
	}
    }

    return wantarray ? sort @files : [sort @files];
}


sub keep_comments{
    my ($self) = @_;
    return $self->{keep_comments} unless wantarray; 
    $self->{keep_comments} = 1;
}

sub no_pair_values {
    my ($self) = @_;
    return not $self->{pair_values} unless wantarray;
    $self->{pair_values} = 0;
}

sub keep_tokens {
    my ($self) = @_;
    return $self->{keep_tokens} unless wantarray;
    $self->{keep_tokens} = 1;
}


# adds files to the files instance var.
sub _set_files_instance_var {
    my ($self, $files) = @_;

    if (defined $files and 'ARRAY' eq ref $files) {

	push @{$self->{files}}, @$files;
	
	$files;
    }
}

# sets sublist to the subs instance var.
sub _set_subs_instance_var {
    my ($self, $subs) = @_;
    if (defined $subs and 'ARRAY' eq ref $subs) {
	$self->{subs} = $subs;  
    }
}


#### Internal helper methods ####


sub _get_perfect_sub_code {
    my ($self) = @_;
    my %cced_sub_code;
    my %perfect_sub_groups = get_perfect_cc_subs($self);
    for my $sub_name (keys %perfect_sub_groups) {    
	my $highest = 0;
	my @matched_files;
	for my $group (@{$perfect_sub_groups{$sub_name}}) {		
	    my $file = @$group[0];				
	    my $code = $self->code($file, $sub_name);		    
	    $code =~ s/./#@$group\ns/;                           #comment specifying files that it came from
	    if (exists $cced_sub_code{join '_', $sub_name, 1}) {  		
		my $new_num = 1 + $highest;
		$code =~ s/sub ([^\s(]+)/sub $1_$new_num/;					
		$cced_sub_code{join '_', $sub_name, $new_num} = $code;
		$highest++;			
	    }		
	    elsif (exists $cced_sub_code{$sub_name}) { 	   	  #second time only
		my $new_code = delete $cced_sub_code{$sub_name}; 		    
		$new_code =~ s/sub ([^\s(]+)/sub $1_1/;		
		$cced_sub_code{join '_', $sub_name, 1} = $new_code;			                       
		redo;
	    } 
	    else {                                                #first time
		$cced_sub_code{$sub_name} = $code;
		$highest = 1;
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

sub _remove_token_keys {
    my ($self) = @_;
    for my $file ($self->files) {
	for my $sub_name (keys %{$self->{data}{$file}}) {
	    delete $self->{data}{$file}{$sub_name}{'tokens'};
	}
    }
}


sub _extract_sub_names {
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
    return \@sub_names;
}


sub _load_files {
    my ($self, @files) = @_;
    foreach my $file (@files) {   
        # Load a document from a file
        my $document = PPI::Document->new($file);	
        # Strip out comments and documentation
	$document->prune('PPI::Token::Pod');
        $document->prune('PPI::Token::Comment')
	    unless $self->keep_comments;

        # Find all the named subroutines
        my $sub_nodes = $document->find(
	    sub { $_[1]->isa('PPI::Statement::Sub') and $_[1]->name } );
	
        next unless $sub_nodes;

        my %sub_names = map {
            $_->name => {
                'code'   => $_->content,
                'tokens' => _token_frequency($_)
            }
        } @$sub_nodes;

        my $file_name = fileparse $file;
        $self->{data}{$file_name} = \%sub_names;
    }

    return $self;
}




sub _cross_correlation {
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
sub _token_frequency {
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














1;






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
