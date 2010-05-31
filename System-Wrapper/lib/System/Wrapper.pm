package System::Wrapper;

use warnings;
use strict;
use overload q{""}    => \&command;
use constant MAX_RAND => 2 ** 32;

use Carp;
use Data::Dumper;
use File::Spec;

our $VERSION = '0.0.1';

sub _fields {
    {   interpreter => undef,
        executable  => undef,
        arguments   => undef,
        input       => undef,
        output      => undef,
        capture     => undef,
        description => undef,
        progress    => undef,
        order       => [qw/interpreter executable arguments input output/],
        path        => [ grep $_, File::Spec->path, q{.} ],
        _destroy =>
            undef,    # internal: safe for additional cleanup in destructor
        _fifo =>
            undef
        ,    # internal: integer unique to current object's fifo pipe if any
        _tmp => int rand
            MAX_RAND,   # internal: integer unique to current object's address
        _arguments =>
            undef,      # internal: stringified form of arguments structure
        _input => undef, # internal: stringified form of input files array ref
        _output =>
            undef
        ,    # internal: stringified form of output file specification hash
    };
}


## CLASS METHODS
sub new {
    my ( $class, %args ) = @_;

    my $self = bless( _fields(), $class );

    while ( my ( $key, $value ) = each %args ) {
        _err( "can't access field '%s' in class '%s'", $key, $class )
            if not exists $self->{$key}
                or $key =~ m/^_/;

        $self->$key($value);
    }
    return $self;
}

sub pipeline {
    _parallel( shift, 'pipe', @_ );
}

sub parallel {
    _parallel( shift, 0, @_ );
}


## OBJECT METHODS
sub interpreter {
    my ( $self, $interpreter ) = @_;

    if ( defined $interpreter ) {
        $self->{interpreter} = $interpreter;
    }

    $self->{_interpreter} = $self->_program_in_path( $self->{interpreter}, 1 );

    return $self->{_interpreter} || $self->{interpreter};
}

sub executable {
    my ( $self, $executable ) = @_;

    if ( defined $executable ) {
        $self->{executable} = $executable;
    }

    $self->{_executable} = $self->_program_in_path( $self->{executable} );

    return $self->{_executable} || $self->{executable};
}

sub arguments {
    my ( $self, $args ) = @_;

    if ( defined $args ) {
        $self->{arguments} = $args;
    }

    $self->{_arguments} = $self->_flatten( $self->{arguments} )
        if $self->{arguments};

    return wantarray ? @{ $self->{arguments} ||= [] } : $self->{_arguments};
}

sub input {
    my ( $self, $input ) = @_;

    if ( defined $input ) {
        _err(
            "type of arg 1 to 'input' must be string or reference to array of strings, not %s",
            ref $input
            )
            unless 'ARRAY' eq ref $input
                or 'SCALAR' eq ref \$input;

        $self->{input} = ref $input ? [$input] : $input;
        $self->{_input} = $self->_flatten( $self->{input} );
    }

    return wantarray ? @{ $self->{input} ||= [] } : $self->{_input};
}

sub output {
    my ( $self, $output ) = @_;

    if ( defined $output ) {
        _err(
            "type of arg 1 to 'output' must be reference to hash with keys 'spec' => 'file', not %s",
            ref $output ? ref $output : ref \$output
        ) unless 'HASH' eq ref $output;

        $self->{output}  = $output;
        $self->{_output} = $self->_flatten( $self->{output} );
    }

    return wantarray ? %{ $self->{output} ||= {} } : $self->{_output};
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
    my ($self, $order) = @_;

    if (defined $order) {
        if ('ARRAY' eq ref $order) {
            if (@$order) {
                _err( "parameter 1 to 'order' must have at most %s elements, not %s",
                      scalar @{$self->{order}}, scalar @$order )
                if scalar @{$self->{order}} < scalar @$order;

                my %spec = map { $_ => 1 } @{$self->{order}};

                for (0 .. @{$self->{order}} - 1) {
                    _err( 
                        "parameter 1 values to 'order' must be any of [%s], not '%s'",
                        join (q{, }, sort keys %spec), $order->[$_] )
                    unless exists $spec{$order->[$_]};
                }
                $self->{order} = $order;
            }
            else {
                _err( "type of parameter 1 to 'order' must be ARRAY, not %s",
                      ref $order || ref \$order );
            }
        }
    }
    return wantarray ? @{$self->{order}} : $self->{order};
}

sub command {
    my ($self) = @_;

    my @command = grep { $_ }
    map { scalar $self->$_ } $self->order;

    return wantarray ? @command : "@command";
}

sub run {
    my ($self) = @_;

    $self->_can_run or return;

    my @command = $self->command;

    #print STDERR q{# }, $self->description, "\n"
    #if $self->description;

    #print STDERR scalar $self->command, "\n";

    if ( $self->output and not -p ( $self->output )[1] ) {
        $command[-1] .= ".part$self->{_tmp}";
    }

    my $stdout;
    eval { $stdout = $self->capture ? qx/"@command"/ : system "@command" };

    $self->_did_run( $@, $? ) or return;

    $self->_rename if $self->output and not -p ( $self->output )[1];

    $self->{_destroy} = 1;

    return $stdout;
}

## INTERNAL METHODS
sub _progress {
    my ($self) = @_;

    return unless $self->progress;

    _err("need input to be set to track progress")
        unless $self->input;

    my $input_size = 0;

    $input_size += -s $_ for $self->input;
}

sub _connect {
    my ( $upstream, $downstream ) = @_;
    my $class = ref $downstream;

    my $named_pipe = File::Spec->catfile(
        File::Spec->tmpdir(),
        qq{$class:} . int rand MAX_RAND
    );

    POSIX::mkfifo( $named_pipe, 0777 )
      or _err( "couldn't create named pipe %s: %s", $named_pipe, $! );

    my %output_spec = $upstream->output;

    _err(
        "can't install fifo %s as output to downstream because there are multiple output specifications (%s):\n%s",
        $named_pipe,
        join( q{, }, map {qq{'$_' }} sort keys %output_spec ),
        $upstream->description || "$upstream"
    ) if keys %output_spec > 1;

    $upstream->output(
        {
            scalar each %output_spec || q{>} => $named_pipe } );
    $downstream->input($named_pipe);

    return $named_pipe;
}

sub _parallel {
    my ( $class, $pipe, $previous, @commands ) = @_;

    return unless @commands;

    eval {
        require POSIX;
        POSIX->import();
        require threads;
    };
    _err( "requires POSIX::mkfifo and threads:\n%s", $@ )
        if $@;

    my @results;
    for my $command ( @commands, 'dummy' ) {
        _err( "type of args to '%s' must be '%s', not '%s'",
            (split /:+/, _this_sub_name(2))[1], $class, ref $command || $command )
            unless ref $command eq $class
                or 'dummy' eq $command;

        if ($pipe) {
            _connect( $previous, $command ) if ref $command;
            $previous->{_fifo} = int rand MAX_RAND;
        }

        # unless (my $pid = fork) {
        #     push @results, $previous->run;
        #     exit;
        # }

        push @results, threads->new( sub { $previous->run } );

        $previous = $command;
    }
    return map { $_->join } @results;
}

sub _rename {
    my ($self) = @_;

    rename $self->_canonical( $self->output . ".part$self->{_tmp}" ),
        $self->_canonical( scalar $self->output )
        or _err(
        "can't rename %s to %s: %s",
        $self->_canonical( $self->output . ".part$self->{_tmp}" ),
        $self->_canonical( scalar $self->output ), $!
        );
}

sub _can_run {
    my ($self) = @_;

    _err("need interpreter or executable to be set to run")
        unless $self->interpreter
            or $self->executable;

    _err( "need interpreter '%s' or executable '%s' to be in path to run",
        $self->interpreter, $self->executable )
        unless $self->{_interpreter}
            or $self->{_executable};

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

sub _program_in_path {
    my ( $self, $program, $is_executable ) = @_;

    return unless $program;

    my ($vol, $dir, $file) = File::Spec->splitpath( $program );

    return $program
	if $dir 
	  and ( -f $program or -l $program )
	    and ( not defined $is_executable or -x $program);

    for ( @{ $self->path } ) {
        my $path = File::Spec->catfile( $_, $program );

        return $path
            if ( -f $path or -l $path )
	      and ( not defined $is_executable or -x $path);
    }

    return q{};
}

sub _canonical {
    my ( $self, $name ) = @_;

    $name =~ s/\s*[<>]\s*//g;
    return File::Spec->canonpath($name);
}

sub _flatten {
    my ( $self, $struct ) = @_;

    return $struct unless ref $struct;

    eval {
	require Storable;
	Storable->import();
    };
    _err( "Storable module required: $@" ) if $@;

    $struct = Storable::dclone($struct);

    my @expanded;
    while ( ref $struct ) {

        if ( 'ARRAY' eq ref $struct ) {
            #$struct = [@$struct];
            last unless @$struct;
            push @expanded, $self->_flatten( shift @$struct );
        }

        elsif ( 'HASH' eq ref $struct ) {
            #$struct = {%$struct};
            last unless %$struct;
            my ( $key, $value ) = each %$struct;
            return unless defined $key and defined $value;
            push @expanded, $key, $self->_flatten($value);
            delete $struct->{$key};
        }

        elsif ( 'SCALAR' eq ref $struct ) {
            #$struct = \${$struct};
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

    #return unless $self->{_destroy};

    if ( $self->{_tmp} ) {
        my ( undef, $out_dir, undef )
            = File::Spec->splitpath( $self->output );
        my $out_files = File::Spec->catfile( $out_dir, "*$self->{_tmp}" );

        unlink glob $out_files;
    }

    if ( $self->{_fifo} ) {
        my $tmp_dir = File::Spec->tmpdir();
        my $tmp_files = File::Spec->catfile( $tmp_dir, "*$self->{_fifo}" );

        unlink glob $tmp_files;
    }
}

## INTERNAL SUBROUTINES
sub _this_sub_name {
    return ( caller( shift || 1 ) )[3];
}

sub _err {
    my ( $spec, @args ) = @_;
    croak sprintf "%s error: $spec",  _this_sub_name(2), @args;
}

1;    # Magic true value required at end of module


package main;
use Carp;
use Data::Dumper;

my $cat = System::Wrapper->new(
    interpreter => 'perl',
    arguments   => [ -pe => q{''} ],
    input       => ['/work/genomes/AT/TAIR_reference.fas.rc'],
    output      => { '>' => 'forward' },
    description =>
        'Concatenate Arabidopsis thaliana reference genome to STDOUT',
);

my $reverse = System::Wrapper->new(
    interpreter => 'perl',
    arguments   => [ -pe => q{'$_ = reverse $_'} ],
    description => 'Reverse input',
);

my $complement = System::Wrapper->new(
    interpreter => 'perl',
    arguments   => [ -pe => q{'tr/ACGT/TGCA/'} ],
    description => 'Reverse input',
    output      => { '>' => 'complement' },
);

my $prog = <<'EOF';
my ($name, $input_size) = (shift, shift);

$input_size = 0 if -f $input_size;

eval {
    require Term::ProgressBar::Simple;
};
die sprintf "%s: requires Term::ProgressBar::Simple:\n%s",
, $@ if $@;

$input_size += -s $_ for @ARGV;

($ENV{COLUMNS}, $ENV{LINES}) = (160, 24);

my $progress = Term::ProgressBar::Simple->new( {
    name  => $name,
    count => $input_size,
    ETA   => "linear",
    fh    => \*STDERR,
} );

use constant BUFFER_SIZE => 4096;
my $buffer               = q{};
my $total_read           = 0;

my $INH = \*STDIN;
while (my $file = shift @ARGV) {

    open $INH, q{<}, $file or die qq{Cant read $file: $!};

    while (my $read = sysread  $INH, $buffer, BUFFER_SIZE) {
        if (syswrite STDOUT, $buffer) {
            $total_read += $read;
            $progress   += $total_read unless int(rand(100000));
        }
    }
    print STDERR "\n";
}
EOF

my $tpb = System::Wrapper->new(
    interpreter => 'perl',
    arguments   => [
        '-e' => qq{'$prog'},
        q{'My Progress Viewer'}, q{121223234}
    ],
    input => ['/work/genomes/AT/TAIR_reference.fas.rc'],
);

my $pv = System::Wrapper->new(
    executable => 'pv',
    arguments  => [ '-s 242427548', '-c -N Test' ],
    input      => ['/work/genomes/AT/TAIR_reference.fas.rc'],
);

croak "failed to run at least one command"
    if grep { $_ }
    System::Wrapper->pipeline( $pv, $cat, $pv, $reverse, $pv, $complement );

__END__


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
    $command->capture(1);
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
880:	final indentation level: 1

Final nesting depth of '('s is 1
The most recent un-matched '(' is on line 409
409: _err(  "%s: $spec", _this_sub_name(2), @args;
         ^
880:	To save a full .LOG file rerun with -g
