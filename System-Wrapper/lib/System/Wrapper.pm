package System::Wrapper;

use warnings;
use strict;
use overload q{""}    => \&stringify;
use constant MAX_RAND => 2 ** 32;

use Carp;
use File::Spec;

our $VERSION = '0.0.1';

## CLASS METHODS
sub new {
    my ( $class, %args ) = @_;
    my $self = bless( { map { $_ => undef }
                        qw/
                              interpreter  executable  arguments  input  output
                              description  path        order
                              progress     spec        capture
                              _interpreter _executable _arguments
                              _input       _output     _destroy
                              _fifo        _tmp        _input_type _output_type
                              debug        verbose
                          /
                      }, $class);

    return $self->_init(%args);
}

## OBJECT METHODS
sub _interpreter : lvalue { shift->{_interpreter} }
sub _executable  : lvalue { shift->{_executable}  }
sub _input       : lvalue { shift->{_input}       }
sub _output      : lvalue { shift->{_output}      }
sub _arguments   : lvalue { shift->{_arguments}   }
sub _fifo        : lvalue { shift->{_fifo}        }
sub _tmp         : lvalue { shift->{_tmp}         }
sub _destroy     : lvalue { shift->{_destroy}     }
sub _input_type  : lvalue { shift->{_input_type}  }
sub _output_type : lvalue { shift->{_output_type} }
sub verbose      : lvalue { shift->{verbose}      }
sub debug        : lvalue { shift->{debug}        }

sub interpreter {
    my ( $self, $interpreter ) = @_;

    if ( defined $interpreter ) {
        $self->{interpreter} = $interpreter;
    }

    $self->_interpreter = $self->_program_in_path( $self->{interpreter}, 1 );

    return $self->_interpreter || $self->{interpreter};
}

sub executable {
    my ( $self, $executable ) = @_;

    if ( defined $executable ) {
        $self->{executable} = $executable;
    }

    $self->_executable = $self->_program_in_path( $self->{executable} );

    return $self->_executable || $self->{executable};
}

sub arguments {
    my ( $self, $args ) = @_;

    if ( defined $args ) {
        $self->{arguments} = $args;
        $self->_arguments  = $self->_flatten( $self->{arguments} )
    }

    return wantarray ? @{ $self->{arguments} ||= [] } : $self->_arguments;
}

sub input {
    my ( $self, $input ) = @_;
    return $self->_io_spec( $input, 'input' );
}

sub output {
    my ( $self, $output ) = @_;
    return $self->_io_spec( $output, 'output' );
}

sub _io_spec {
    my ($self, $spec, $io) = @_;

    my $private_method = "_$io";

    if ( defined $spec ) {
        
        my $io_type = $private_method . '_type';
        $self->$io_type        = ref $spec || 'SCALAR';
        $self->{$io}           = $spec;
        $self->$private_method = $self->_flatten( $self->{$io} );
    }

    return wantarray ? $self->_deref( $self->{$io} ) : $self->$private_method;
}

sub path {
    my ( $self, $path ) = @_;

    if ( defined $path ) {
        $path = [$path] if 'SCALAR' eq ref \$path;

        for (@$path) {
            _err(
                "type of arg 1 to 'path' must be directory or reference to array of directories (not non-directory '%s')",
                $_
            ) unless -d $_;
        }

        $self->{path} = $path;
    }
    return wantarray ? @{ $self->{path} ||= [] } : $self->{path};
}

sub capture {
    my ( $self, $capture ) = @_;

    if ($capture) {
        $self->{capture} = $capture;
    }
    return $self->{capture};
}

sub description {
    my ( $self, $description ) = @_;

    if ($description) {
        $self->{description} = $description;
    }
    return $self->{description};
}

sub progress {
    my ( $self, $progress ) = @_;

    if ($progress) {
        $self->{progress} = $progress;
    }
    return $self->{progress};
}

sub order {
    my ( $self, $order ) = @_;

    if ( defined $order ) {
        if ( 'ARRAY' eq ref $order ) {
            if (@$order) {
                _err(
                    "parameter 1 to 'order' must have at most %s elements, not %s",
                    scalar @{ $self->{order} },
                    scalar @$order
                ) if scalar @{ $self->{order} } < scalar @$order;

                my %spec = map { $_ => 1 } @{ $self->{order} };

                for ( 0 .. @{ $self->{order} } - 1 ) {
                    _err(
                        "parameter 1 values to 'order' must be any of [%s], not '%s'",
                        join( q{, }, sort keys %spec ),
                        $order->[$_]
                    ) unless exists $spec{ $order->[$_] };
                }
                $self->{order} = $order;
            }
            else {
            }
        }
        else {
            _err( "type of parameter 1 to 'order' must be ARRAY, not %s",
                ref $order || ref \$order );
        }
    }
    return wantarray ? @{ $self->{order} } : $self->{order};
}

sub command {
    my ($self) = @_;

    my @command = grep {$_}
        map { scalar $self->$_ } $self->order;

    return wantarray ? @command : "@command";
}

sub stringify {
    my ($self) = @_;

    my $out = $self->description ? q{# } . $self->description . qq{\n}
                                 : q{}
            . $self->command                                  . qq{\n};

    return $out;
}

sub run {
    my ($self) = @_;
    
    if ($self->debug) {
        $self->_rename_to_tmp;
        print STDERR $self->stringify;
        $self->_rename_from_tmp;
        return;
    }

    if ($self->verbose) {
        print STDERR $self->stringify;
    }

    $self->_can_run;

    $self->_rename_to_tmp;

    my @command = $self->command;

    my $stdout;
    eval { $stdout = $self->capture ? qx/"@command"/ : system "@command" };

    $self->_did_run( $@, $? ) or return;

    $self->_rename_from_tmp;

    $self->_destroy = 1;

    return $stdout;
}

## INTERNAL METHODS
sub _init {
    my ( $self, %args ) = @_;

    $self->{order} = [qw/interpreter executable arguments input output/];
    $self->{path}  = [ grep {$_} File::Spec->path, q{.} ];
    $self->{_tmp}  = int rand MAX_RAND;

    while ( my ( $key, $value ) = each %args ) {

        _err( "Can't access '%s' field in class '%s'",
              $key, ref $self )
        if !exists $self->{$key}
        or $key =~ m/^_/;

        $self->$key($value);
    }

    return $self;
}

sub _progress {
    my ($self) = @_;

    return unless $self->progress;

    _err("need input to be set to track progress")
        unless $self->input;

    my $input_size = 0;

    $input_size += -s $_ for $self->input;
}

sub _tmp_name {
    my ($self, $name) = @_;
    return $name . '.part' . $self->_tmp;
}

sub _rename_to_tmp {
    my ($self) = @_;

    return unless $self->output;

    if ('HASH' eq $self->_output_type) {
        my %output = $self->output;
      NAME:
        for my $spec ( keys %output ) {

            use Data::Dumper;
            print STDERR Dumper \%output if not defined $output{$spec};

            next NAME if -p $output{$spec};
            my $tmp = $self->_tmp_name($output{$spec});

            $output{$spec} = $tmp;
        }
        $self->output(\%output);
    }
    elsif ('ARRAY' eq $self->_output_type) {
        my @output = $self->output;
      NAME:
        for my $name (@output) {
            next NAME if -p $name;
            my $tmp = $self->_tmp_name($name);
            
            $name = $tmp;
        }
        $self->output(\@output);
    }
    elsif ('SCALAR' eq $self->_output_type) {
        my ($name) = $self->output;
        my $tmp = $self->_tmp_name($name);
        $self->output($tmp);
    }
    else {
        _err( 
            "%s is not an accepted output type. It should be impossible to get this error",
            $self->_output_type
        );
    }
}

sub _rename_from_tmp {
    my ($self) = @_;

    return unless $self->output;

    my (@name, @tmp);

    if ('HASH' eq $self->_output_type) {
        my %output = $self->output;
      NAME:
        for my $spec ( keys %output ) {
            push @tmp,  $output{$spec};
            push @name, $output{$spec};
            $name[-1] =~ s/\.part.+//;
        }
    }
    elsif ('ARRAY' eq $self->_output_type) {
        my @output = $self->output;
      NAME:
        for my $tmp (@output) {
            push @tmp, $tmp;
            push @name, $tmp;
            $name[-1] =~ s/\.part.+//;
        }
    }
    elsif ('SCALAR' eq $self->_output_type) {
        push @tmp, $self->output;
        push @name, $tmp[-1];
        $name[-1] =~ s/\.part.+//;
    }
    else {
        _err( 
            "%s is not an accepted output type. It should be impossible to get this error",
            $self->_output_type
        );
    }

  NAME:
    for (0 .. @name - 1) {
        my ($name, $tmp) = ($name[$_], $tmp[$_]);
        next NAME if -p $name;

        my $mv = System::Wrapper->new(
            interpreter => 'perl',
            arguments   => qq{-e 'rename "$tmp", "$name" or die "\$!"'},
        );
        $mv->verbose = 1 if $self->verbose;
        $mv->debug   = 1 if $self->debug;
        $mv->run;
    }
}

sub _can_run {
    my ($self) = @_;

    if ( my $nonexistent_file = $self->_inputs_not_available ) {
        _err( "%s %s", $nonexistent_file,
            -e $nonexistent_file
            ? 'is not readable'
            : 'does not exist' );
    }

    _err("need interpreter or executable to be set to run")
        unless $self->interpreter
            or $self->executable;

    _err( "need interpreter '%s' or executable '%s' to be in path to run",
        $self->interpreter, $self->executable )
        unless $self->_interpreter
            or $self->_executable;

    return 1;
}

sub _did_run {
    my ( $self, $eval_error, $child_error ) = @_;

    _err( "failed to run command:\n%s\n%s",
        scalar $self->command, $eval_error )
        if $eval_error;

    _err(
        "failed to run command (%s):\n%s",
        ( $child_error & 255 )
        ? 'signal ' . ( $child_error & 255 )
        : 'exit ' . ( $child_error >> 8 ),
        scalar $self->command
    ) if $child_error;

    return 1;
}

sub _inputs_not_available {
    my ($self) = @_;

    for ( $self->input ) {
        return $_ unless -r $_;
    }

    return;
}

sub _program_in_path {
    my ( $self, $program, $is_executable ) = @_;

    return unless $program;

    my ( $vol, $dir, $file ) = File::Spec->splitpath($program);

    return $program
        if $dir
            and ( -f $program or -l $program )
            and ( not defined $is_executable or -x $program );

    for ( @{ $self->path } ) {
        my $path = File::Spec->catfile( $_, $program );

        return $path
            if ( -f $path or -l $path )
            and ( not defined $is_executable or -x $path );
    }

    return q{};
}

sub _deref {
    my ($self, $struct) = @_;

    return unless $struct;
    return $struct unless ref $struct;
    return $$struct if 'SCALAR' eq ref $struct;
    return @$struct if 'ARRAY'  eq ref $struct;
    return %$struct if 'HASH'   eq ref $struct;
}

sub _flatten {
    my ( $self, $struct ) = @_;

    return $struct unless ref $struct;

    eval {
        require Storable;
        Storable->import();
    };
    _err("Storable module required: $@") if $@;

    $struct = Storable::dclone($struct);

    my @expanded;
    while ( ref $struct ) {

        if ( 'ARRAY' eq ref $struct ) {
            last unless @$struct;
            push @expanded, $self->_flatten( shift @$struct );
        }

        elsif ( 'HASH' eq ref $struct ) {
            last unless %$struct;
            my ( $key, $value ) = each %$struct;
            return unless defined $key and defined $value;
            push @expanded, $key, $self->_flatten($value);
            delete $struct->{$key};
        }

        elsif ( 'SCALAR' eq ref $struct ) {
            last unless $$struct;
            push @expanded, $self->_flatten($$struct);
            $struct = undef;
        }
        else {
            _err(
                "type of arg 1 to '_flatten' must be a scalar, hash or array reference (not '%s')",
                ref $struct
            );
        }
    }
    @expanded = grep $_, @expanded;
    return wantarray ? @expanded : "@expanded";
}

sub DESTROY {
    my ($self) = @_;

    return unless $self->_destroy;

    if ( $self->_tmp ) {
        my ( undef, $out_dir, undef )
            = File::Spec->splitpath( $self->output );
        my $out_files = File::Spec->catfile( $out_dir, q{*} . $self->_tmp );

        unlink glob $out_files;
    }

    if ( $self->_fifo ) {
        unlink glob $self->_fifo;
    }
}

## INTERNAL SUBROUTINES
sub _this_sub_name {
    return ( caller( shift || 1 ) )[3];
}

sub _err {
    my ( $spec, @args ) = @_;
    croak sprintf "%s error: $spec", _this_sub_name(2), @args;
}


1;    # Magic true value required at end of module


=head1 NAME

System::Wrapper - Class-wrapped system calls and qx operator

=head1 VERSION

This document describes System::Wrapper version 0.0.1


=head1 SYNOPSIS

    use System::Wrapper;

    my $command = System::Wrapper->new();

    $command->interpreter( 'perl');
    $command->executable( 'program.pl');
    $command->arguments( [ 'first', {second_a => 2, second_b => 2}, {third => [1,2,3]} ] );
    $command->input( \@ARGV );
    $command->output( { '--output' => 'file'}, q{>} => 'file2' );
    $command->path( [$command->path, q{.}] );
    $command->capture = 1;
    $command->verbose = 1;
    print $command->command;
    $command->run;


=head1 DESCRIPTION

 This module wraps perl's C<system> call and c<qx> operator in an object-oriented
 interface. It provides utility methods for accomplishing things that are not very
 simple in C<system> and C<qx>. This includes in-situ I/O and call success via
 temporary filenames, C<system> call progress estimation, finding whether the
 executable and-or interpreter are on the path, validating filenames, cross-platform
 output operators and argument type specification.

 This module can be used as a generic wrapper around C<system> and C<qx>, or as
 a base class for building interfaces to utilities not available to C<perl> itself.

=head1 INTERFACE 

=head2 CLASS METHODS

=over

=item new(%args)

    my %args = (
        interpreter => undef, # optional: string
        executable  => undef, # required: string
        arguments   => undef, # optional: any nested structure of hashes,
                              # arrays or scalar references
        input       => undef, # optional: scalar or array reference
        output      => undef, # optional: hash reference of form { spec => file }
                              # eg:   { '>' => 'out' } or { '--output' => 'out' }
        capture     => undef, # optional: return stdout, instead of exit code,
                              # via $self->run
        path        => [ grep $_, File::Spec->path, q{.} ]
                              # required: path of directories on which to look for
                              # interpreter and executable programs
    );

    my $command = System::Wrapper->new(%args);

=back


=head2 SELECTOR METHODS

new
interpreter
executable
arguments
input
output
path
capture
command
run
_program_in_path
_canonical
_flatten
DESTROY

=over

=item server_uri()

=item server_uri($uri)

Default C<$uri>: L<http://api.wordnik.com/api-v3>


=item api_key()

=item api_key($key)

Required C<$key>: Your API key, which can be requested at L<http://api.wordnik.com/signup/>.


=item version()

=item version($version)

Default C<$version>: I<3>. Only API version 3 (the latest) is currently supported.


=item format()

=item format($format)

Default C<$format>: I<json>. Other accepted formats are I<xml> and I<perl>.


=item cache()

=item cache($cache)

Default C<$cache>: I<10>. Number of requests to cache. Deletes the oldest request if cache fills up.


=item debug()

=item debug($debug)

Default C<$debug>: I<0>. Don't sent GET requests to Wordnik. Return the actual request as a string.

=back


=head2 OBJECT METHODS

=over

=item word($word, %args)

This returns the word you requested, assuming it is found in our corpus.
See L<http://docs.wordnik.com/api/methods#words>.

C<$word> is the word to look up. C<%args> accepts:

Default C<useSuggest>: I<false>. Return an array of suggestions, if available.

Default C<literal>: I<true>. Return non-literal matches.

If the suggester is enabled, you can tell it to return the best match with C<useSuggest=true> and C<literal=false>.


=item phrases($word, %args)

You can fetch interesting bi-gram phrases containing a word.
The "mi" and "wlmi" elements refer to "mutual information" 
and "weighted mutual information" and will be explained in detail via future blog post.
See L<http://docs.wordnik.com/api/methods#phrases>.

C<$word> is the word to look up. C<%args> accepts:

Default C<count>: I<5>. Specify the number of results returned.


=item definitions($word, %args)

Definitions for words are available from Wordnik's keying of the Century Dictionary and parse of the Webster GCIDE.
The Dictionary Model XSD is available in L<http://github.com/wordnik/api-examples/blob/master/docs/dictionary.xsd> in GitHub.
See L<http://docs.wordnik.com/api/methods#definitions>.

C<$word> is the word to look up. C<%args> accepts:

Default C<count>: I<5>. Specify the number of results returned.

Default C<partOfSpeech>: I<empty>. Specify one or many part-of-speech types for which to return definitions. Pass multiple types as an array reference.

The available partOfSpeech values are:

    [noun, verb, adjective, adverb, idiom, article, abbreviation, preposition, prefix, interjection, suffix]


=item examples($word)

You can retrieve 5 example sentences for a words in Wordnik's alpha corpus. Each example contains the source document and a source URL, if it exists.
See L<http://docs.wordnik.com/api/methods#examples>.

C<$word> is the word to look up.


=item related($word, %args)

You retrieve related words for a particular word.
See L<http://docs.wordnik.com/api/methods#relateds>.

C<$word> is the word to look up. C<%args> accepts:

Default C<type>: I<empty>. Return one or many relationship types. Pass multiple types as an array reference.

The available type values are:

    [synonym, antonym, form, equivalent, hyponym, variant]


=item frequency($word)

You can see how common particular words occur in Wordnik's alpha corpus, ordered by year.
See L<http://docs.wordnik.com/api/methods#freq>.

C<$word> is the word to look up.


=item punctuationFactor($word)

You can see how common particular words are used with punctuation.
See L<http://docs.wordnik.com/api/methods#punc>.

C<$word> is the word to look up.


=item suggest($word, %args)

The autocomplete service gives you the opportunity to take a word fragment (start of a word) and show what other words start with the same letters.
The results are based on corpus frequency, not static word lists, so you have access to more dynamic words in the language.
See L<http://docs.wordnik.com/api/methods#auto>.

C<$word> is the word to look up. C<%args> accepts:

Default C<count>: I<5>. Specify the number of results returned.

Default C<startAt>: I<0>. You can also specify the starting index for the results returned. This allows you to paginate through the matching values.


=item wordoftheday

You can fetch Wordnik's word-of-the day which contains definitions and example sentences.
See L<http://docs.wordnik.com/api/methods#wotd>.


=item randomWord(%args)

You can fetch a random word from the Alpha Corpus.
See L<http://docs.wordnik.com/api/methods#random>.

C<%args> accepts:

Default C<hasDictionaryDef>: I<true>. You can ask the API to return only words where there is a definition available.

=back


=head1 INSTALLATION

To install this module type the following:

   perl Build.PL
   Build
   Build test
   Build install

or

   perl Makefile.PL
   make
   make test
   make install


=head1 DIAGNOSTICS

=over

=item C<< "Can't access '$key' field in class $class" >>

Private or inexistent member variable.

=item C<< "Invalid argument key or value: '$type'" >>

Inexistent query parameter, or wrong value passed to existing parameter.

=item C<< "Parameter 'partOfSpeech' requires a reference to an array" >>

partOfSpeech => [qw/.../].

=item C<< "Parameter 'type' requires a reference to an array" >>

type => [qw/.../].

=item C<< "The operation you requested requires JSON to be installed" >>

perl -MCPAN -e 'install JSON'.

=item C<< "Unsupported api format: '$format'" >>

Supported formats are 'perl', 'json', 'xml'.

=item C<< "Unsupported api version: '$version'" >>

The only API version supported by this module is 3.

=back


=head1 CONFIGURATION AND ENVIRONMENT

System::Wrapper requires no configuration files or environment variables.


=head1 DEPENDENCIES

This module requires the core modules C<Test::More>, C<version> and C<Carp>, and C<LWP::UserAgent> from C<CPAN>.
Additionally, it recommends-requires C<JSON> from C<CPAN> for getting data in Perl data structures.


=head1 INCOMPATIBILITIES

None reported.


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Response headers are not checked for 404s, etc. Likewise, response data is not post-processed in any way, other
than optionally being parsed from C<JSON> to C<Perl> data structures. Data::Dumper should be of help there.

Please report any bugs or feature requests to
C<bug-www-wordnik-api@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.


=head1 TODO

=over

=item Error checking

Implement basic HTTP error checking on response headers.

=item Post-processing

Add filtering methods on response data.

=back


=head1 AUTHOR

Pedro Silva  C<< <pedros@berkeley.edu> >>


=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010, Pedro Silva C<< <pedros@berkeley.edu> >>. All rights reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see L<http://www.gnu.org/licenses/>.
