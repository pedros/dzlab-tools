package DZLab::Tools::GFFStore;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use DBI;
use DZLab::Tools::GFF qw/parse_gff_arrayref/;

my @default_cols     = qw/seqname source feature start   end     score strand frame   attribute/;
my @default_coltypes = qw/text    text   text    numeric numeric real  text   numeric text/;

sub new {
    my $class = shift;
    my $opt = shift;
    my $self = {};

    $self->{dbname}      = $opt->{dbname}     || ':memory:';
    $self->{attributes}  = [];
    $self->{indices}     = $opt->{indices}    || [];
    $self->{debug}       = $opt->{debug}      || 0;
    $self->{verbose}     = $opt->{verbose}    || 0;
    $self->{columns}     = [@default_cols];
    $self->{columntypes} = [@default_coltypes];
    $self->{counter}     = 10000;


    while (my ($col,$coltype) = each %{$opt->{attributes}}) {
        push @{$self->{attributes}}, $col;
        push @{$self->{columns}}, $col;
        push @{$self->{columntypes}}, $coltype;
    }

    $self->{numcol}      = scalar @{$self->{columns}};

    my $blessed = bless $self, $class;

    #$blessed->slurp();

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

sub slurp{
    my $self = shift;
    my $opt = shift;

    my $filename = $opt->{filename}  || undef;
    my $handle   = $opt->{handle}    || undef;
    my $append   = $opt->{append} || undef;
    unless ($filename xor $handle) {croak "Need handle or a filename"};

    my $dbname = $self->{dbname};
    unlink $dbname if (!$self->{append} && $dbname ne ':memory:');

    # create database
    my $dbh = $self->{dbh} = DBI->connect("dbi:SQLite:dbname=$dbname","","",
        {RaiseError => 1, AutoCommit => 1});
    $dbh->do("PRAGMA automatic_index = OFF");
    $dbh->do("PRAGMA journal_mode = OFF");
    $dbh->do($self->create_table_statement());

    #say $self->insert_statement();
    my $insert_sth = $dbh->prepare($self->insert_statement());

    $dbh->{AutoCommit} = 0;

    my $fh;
    if ($filename){
        say "opening $filename" if $self->{verbose};
        open $fh, '<', $filename or croak "can't open $filename";
    } else{
        $fh = $handle;
    }
    my $counter = 0;
    while (my $line = <$fh>){
        chomp $line;
        my $parsed = parse_gff_arrayref($line,@{$self->{attributes}});
        next unless $parsed;

        $insert_sth->execute(@$parsed);
        if ($counter++ % $self->{counter} == 0){
            say "Read " . ($counter-1) if $self->{verbose};
            $dbh->commit;
        }
    }
    $dbh->commit;
    $dbh->{AutoCommit} = 1;
    if ($filename){
        close $fh;
    }
}

sub create_indices{
    my $self = shift;
    my $dbh  = $self->{dbh};
    say "creating indices (if any)" if $self->{verbose};
    for my $index_statement ($self->create_index_statements()){
        say $index_statement if $self->{verbose};
        $dbh->do($index_statement);
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

=head2 select_iter

Run raw select statement against db, return an arrayref iterator

=cut

sub select_iter{
    croak("not yet implemented");
}

=head2 select

Run raw select statement against db, return an arrayref of arrayrefs

=cut

sub select{
    croak("not yet implemented");
}

=head2 make_iterator {column1 => value1, column2 => value2, ...}

return an iterator which, for every call, returns a hashref of a row matching the given 
equality constraints.

=cut

sub make_iterator{
    my $self = shift;
    my $constraints = shift;
    my @values;

    my @where;
    while (my ($col,$val) = each %$constraints) {
        if (defined($val)){
            push @where, "$col = ?";
            # keep track of @values in same order as added to where
            push @values, $val;
        } else{
            push @where, "$col is null";
            # $val is not pushed here since does not require a placeholder
        }
    }
    my $where_clause = @where ? ("where " . join " and ", @where) : "";

    my $select_stmt = "select * from gff $where_clause";
    say $select_stmt if $self->{debug};

    my $dbh = $self->{dbh};
    my $sth = $dbh->prepare($select_stmt);
    $sth->execute(@values);
    return sub{
        return $sth->fetchrow_hashref();
    };
}

=head2 query {column1 => value1, column2 => value2, ...}

like make_iterator, except slurps up entire row set and returns arrayref

=cut

sub query{
    my $self = shift;
    my $constraints = shift;
    my $iter = $self->make_iterator($constraints);
    my @accum;
    while (my $row = $iter->()){
        push @accum,$row;
    }
    return \@accum;
}

=head2 exists {column1 => value1, column2 => value2, ...}

=cut

sub exists{
    my $self = shift;
    my $constraints = shift;
    my $iter = $self->make_iterator($constraints);
    my @accum;
    if (my $row = $iter->()){
        return 1;
    }
    return 0;
}

=head2 make_iterator_overlappers [[$start1, $end1], [start2, $end2], ...]

returns and iterator which, on every call, returns a element overlapping with the ranges.
may return same thing twice....

=cut

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

=head2 $gffstore->sequences

return distinct elements from the seqname column

=cut

sub sequences{
    my $self = shift;
    my $dbh = $self->{dbh};
    my $results = $dbh->selectall_arrayref("select distinct seqname from gff");

    return map { $_->[0] } @$results;
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
    $self->{dbh}->disconnect;
}

1;

