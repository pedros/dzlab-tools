package DZLab::Tools::GFFStore;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use DBI;

my @default_cols     = qw/seqname source feature start   end     score strand frame   attribute/;
my @default_coltypes = qw/text    text   text    numeric numeric real  text   numeric text/;

sub new {
    my $class = shift;
    my $opt = shift;
    my $self = {};

    $self->{attributes}  = [];
    $self->{filename}    = $opt->{filename}   || undef;
    $self->{handle}      = $opt->{handle}   || undef;
    $self->{indices}     = $opt->{indices}    || [];
    $self->{debug}       = $opt->{debug}      || 0;
    $self->{verbose}     = $opt->{verbose}    || 0;
    $self->{columns}     = [@default_cols];
    $self->{columntypes} = [@default_coltypes];
    $self->{counter}     = 10000;

    unless ($self->{handle} xor $self->{filename}) {croak "Need handle or a filename"};

    while (my ($col,$coltype) = each %{$opt->{attributes}}) {
        push @{$self->{attributes}}, $col;
        push @{$self->{columns}}, $col;
        push @{$self->{columntypes}}, $coltype;
    }

    $self->{numcol}      = scalar @{$self->{columns}};
    #say $self->{numcol};

    my $blessed = bless $self, $class;

    $blessed->slurp();

    return $blessed;
}
sub insert_statement{
    my $self = shift;
    my $colcomma = join ",", @{$self->{columns}};
    my $placeholders = join (q{,}, map {'?'} (0 .. $self->{numcol} - 1));
    return "insert into gff ($colcomma) values ($placeholders)";
}
sub create_table_statement{
    my $self = shift;
    return "create table gff (" . 
    join(',',
        map { $self->{columns}[$_] . " " .  $self->{columntypes}[$_]} (0 .. $self->{numcol}-1)
    ) 
    .  ")";
}

sub create_index_statements{
    my $self = shift;
    
    return map { 
        my @cols = @$_;
        my $colcomma   = join q{,}, @cols;
        my $index_name = join q{_}, @cols;
        return "create index $index_name on gff ($colcomma)";
    } @{$self->{indices}};
}

sub _get_attr{
    my $attr = shift;
    my %accum = map {$_ => undef} @_;

    my $key;
    while($attr =~ m/
            (?:
            \s*?([^=;]+)\s*?  # key, which may not be there
            =
            )?
            \s*?([^=;]+)\s*?  # value
            /xmsg){
        $key = $1 ? $1 : 'Note';
        if (exists $accum{$key}){
            $accum{$key} = $2;
        }
    }
    return @accum{@_};
}

sub slurp{
    my $self = shift;

    my $dbname = $self->{debug} ? 'debug.db' : ':memory:';
    unlink 'debug.db' if $self->{debug};

    # create database
    my $dbh = $self->{dbh} = DBI->connect("dbi:SQLite:dbname=$dbname","","");
    $dbh->{RaiseError} = 1;
    $dbh->do("PRAGMA automatic_index = OFF");
    $dbh->do("PRAGMA journal_mode = OFF");
    $dbh->do($self->create_table_statement());

    #say $self->insert_statement();
    my $insert_sth = $dbh->prepare($self->insert_statement());

    $dbh->{AutoCommit} = 0;

    my $fh;
    if ($self->{filename}){
        say "opening $self->{filename}" if $self->{verbose};
        open $fh, '<', $self->{filename} or croak "can't open $self->{filename}";
    } else{
        $fh = $self->{handle};
    }
    my $counter = 0;
    while (my $line = <$fh>){
        chomp $line;
        $line =~ s/[\n\r]//g;
        my @arr = split /\t/, $line;
        next unless @arr == 9;
        my @attrvals = _get_attr($arr[8],@{$self->{attributes}});

        $insert_sth->execute(@arr,@attrvals) ;
        if ($counter++ % $self->{counter} == 0){
            say "Read " . ($counter-1) if $self->{verbose};
            $dbh->commit;
        }
    }
    $dbh->commit;
    $dbh->{AutoCommit} = 1;

    say "creating indices (if any)" if $self->{verbose};
    for my $index_statement ($self->create_index_statements()){
        say $index_statement;
        $dbh->do($index_statement);
    }
    
    if ($self->{filename}){
        close $fh;
    }
}

##########################################################
# Accessors

sub count{
    my $self = shift;
    my $dbh = $self->{dbh};
    my @r = $dbh->selectrow_array("select count(*) from gff");
    return $r[0];
}

sub make_iterator_arrayref{
    my $self = shift;
    my $dbh = $self->{dbh};
    my $select = $dbh->prepare("select * from gff");
    $select->execute();
    return sub{
        return $select->fetchrow_arrayref();
    };
}

sub make_iterator{
    my $self = shift;
    my $dbh = $self->{dbh};
    my $select = $dbh->prepare("select * from gff");
    $select->execute();
    return sub{
        return $select->fetchrow_hashref();
    };
}

sub make_iterator_overlappers{
    my $self = shift;
    my $ranges = shift;

    my $dbh = $self->{dbh};

    my @iterators; 
    for my $range (@$ranges){
        push @iterators, $dbh->prepare("select * from gff where start <= $range->[1] and end >= $range->[0]");
    }

    return sub{ #rewrite this
        while (@iterators){
            my $first = $iterators[0];
            $first->execute() unless $first->{Executed};
            my $row = $first->fetchrow_hashref();
            if ($row){
                return $row;
            } else{
                #done with this one
                shift @iterators;
            }
        }
        return;
    };
}

sub dump{
    my $self = shift;
    my $dbh = $self->{dbh};
    my $select = $dbh->prepare("select * from gff");
    $select->execute();
    while (my $row = $select->fetchrow_arrayref()){
        say Dumper $row;
    }
}

sub DESTROY{
    my $self = shift;
    close $self->{dbh};
}

1;

