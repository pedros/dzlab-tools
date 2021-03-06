package DZLab::Tools::GFFStore;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use feature 'state';
use Carp;
use DBI;
use Regexp::Common;
use File::Temp qw/tempfile unlink0/;
use DZLab::Tools::GFF qw/parse_gff_arrayref gff_to_string/;

my @default_cols     = qw/seqname source feature start   end     score strand frame   attribute/;
my @default_coltypes = qw/text    text   text    numeric numeric real  text   numeric text/;


=head2 new

 my $gffstore = DZLab::Tools::GFFStore->new({
    dbname     => 'filename',
    indicies   => [['id'], ['seqname']],
    attributes => {id => 'text', c => 'numeric'},
    });


=head3 options 

By default, temporary on-disk database used. As an alternative, 'memory' for in-memory db, or 'dbname' for a specific
file.

 dbname     - filename of output sqlite database
 memory     - in memory database
 tmp        - use temporary (default)

 overwrite  - for dbname only. whether to delete file first (default = 1)

 indices    - arrayref of index, where each index is an arrayref of column names
 attributes - hashref of attribute field name => field type ('numeric', 'text')

=cut
sub new {
    my $class = shift;
    my $opt = shift;
    #my $self = {};
    #my $blessed = bless $self, $class;
    my $self = bless {}, $class;

    # the three possible sqlite storage options. default to temp file since the 
    # operations seem to be cpu bound more than IO bound. may change after benchmarking.
    $self->{dbname}      = $opt->{dbname}     || 0;
    $self->{memory}      = $opt->{memory}     || 0;
    $self->{tmp}         = $opt->{tmp}        || 0;
    if (!$self->{dbname} && ! $self->{memory} && ! $self->{tmp}){
        $self->{tmp} = 1;
    }

    # only for dbname option: if existing db should be overwritten.
    # default to yes since can't imagine reusing db's (yet)
    $self->{overwrite}   = $opt->{overwrite}     || 1;

    # attributes and indices
    $self->{attributes}  = []; # not an error, this should be array. see line 67
    $self->{indices}     = $opt->{indices}    || [];
    $self->{columns}     = [@default_cols];
    $self->{columntypes} = [@default_coltypes];

    # other 
    $self->{debug}       = $opt->{debug}      || 0;
    $self->{verbose}     = $opt->{verbose}    || 0;
    $self->{counter}     = $opt->{commitsize} || 10000;

    croak "dbname, memory, temp options are mutually exclusive"
    unless ( $self->{dbname} xor $self->{memory} xor $self->{tmp});

    # use the dbname attribute even for memory and tmp
    if ($self->{dbname} && -e $self->{dbname} && $self->{overwrite}){ # overwrite?
        unlink $self->{dbname};
    }
    elsif ($self->{memory}) {
        $self->{dbname} = ':memory:';
    } 
    elsif ($self->{tmp}) {
        ($self->{tmpfh}, $self->{dbname}) = tempfile();
        if($self->{debug}){
            say STDERR "tmpfile = $self->{dbname}";
        }
    }

    # go through attributes and add them to @columns and @columntypes. 
    # keeping their order is important b/c when inserting via param bind, 
    # doesn't accepting hashes
    if (ref $opt->{attributes} eq 'HASH'){ # but only if we were actually given attr hash
        while (my ($col,$coltype) = each %{$opt->{attributes}}) {
            push @{$self->{attributes}}, $col;
            push @{$self->{columns}}, $col;
            push @{$self->{columntypes}}, $coltype;
        }
    }
    $self->{numcol}      = scalar @{$self->{columns}};

    my $dbname = $self->{dbname};

    # create database
    my $dbh = $self->{dbh} = DBI->connect("dbi:SQLite:dbname=$dbname","","",
        {RaiseError => 1, AutoCommit => 1});
    $dbh->do("PRAGMA automatic_index = OFF");
    $dbh->do("PRAGMA journal_mode = OFF");
    $dbh->do("PRAGMA cache_size = 80000");
        
    $dbh->do($self->create_table_statement());

    return $self;
}

sub insert_statement{
    my $self = shift;
    my $colcomma = join ",", @{$self->{columns}};
    my $placeholders = join (q{,}, map {'?'} (0 .. $self->{numcol} - 1));
    return "insert into gff ($colcomma) values ($placeholders)";
}

sub create_table_statement{
    my $self = shift;
    return "create table if not exists gff (_id integer primary key autoincrement, " . 
    join(',',
        map { $self->{columns}[$_] . " " .  $self->{columntypes}[$_]} (0 .. $self->{numcol}-1)
    ) 
    .  ")";
}

sub create_rtree_statement{
    my $self = shift;
    return "create virtual table range using rtree (_id integer primary key, start, end)";
}
sub insert_rtree_statement{
    return "insert into range (_id, start, end) values (?,?,?)";
}

sub create_index_statements{
    my $self = shift;
    
    return (map { 
        my @cols = @$_;
        my $colcomma   = join q{,}, @cols;
        my $index_name = join q{_}, @cols;
        "create index $index_name on gff ($colcomma)";
    } @{$self->{indices}});
}

=head2 slurp

 $gffstore->slurp({handle => \*HANDLE});
 $gffstore->slurp({filename => "genes.gff"});

Read a filehandle or a file into store.

=cut
sub slurp{
    my $self = shift;
    my $opt = shift;

    my $filename = $opt->{filename}  || undef;
    my $handle   = $opt->{handle}    || undef;
    unless ($filename xor $handle) {croak "Need handle or a filename"};

    my $dbh = $self->{dbh};

    say STDERR $self->insert_statement() if $self->{debug};
    my $insert_sth = $dbh->prepare($self->insert_statement());

    $dbh->{AutoCommit} = 0;

    my $fh;
    if ($filename){
        say STDERR "opening $filename" if $self->{verbose};
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
            say STDERR "Reading "  . ($filename ? $filename : q{}) . ' ' . ($counter-1) if $self->{verbose};
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
    say STDERR "creating indices (if any)" if $self->{verbose};

    for my $index_statement ($self->create_index_statements()){
        say STDERR $index_statement if $self->{verbose};
        $dbh->do($index_statement);
    }
}
    

sub rtree_supported{
    return scalar eval{
        my $dbh = DBI->connect("dbi:SQLite:dbname=:memory:","","", {RaiseError => 1});
        $dbh->do("create virtual table x using rtree(x,y,z)");
        1;
    } || 0;
}


##########################################################
# Accessors

=head2 count

 $gffstore->count();

Return the number of rows in store.

=cut
sub count{
    my $self = shift;
    my $dbh = $self->{dbh};
    my @r = $dbh->selectrow_array("select count(*) from gff");
    return $r[0];
}

=head2 select_iter

 my $iter = $gffstore->select_iter("select seqname from gff");

 while (defined(my $row_href = $iter->())){
    ...
 }

Run raw select statement against db, return an hashref iterator

=cut
sub select_iter{
    my $self = shift;
    my $stmt = shift;
    my $dbh = $self->{dbh};
    my $sth = $dbh->prepare($stmt);
    $sth->execute(@_) or die "can't execute statement";
    
    return sub {
        return $sth->fetchrow_hashref();
    };
}


=head2 select_aggregate $feature, $by, @on

Run aggregating select over rows of type $feature grouped by $by on cols @on

Return an arrayref iterator of [[aggregate-cols], group-by-col]

=cut
sub select_aggregate {
    my ($self, $feature, $by, @on) = @_;
    my $on = join ',', @on;

    require DZLab::Tools::ArrayAggregator;
    $self->{dbh}->sqlite_create_aggregate( 'aggregate', -1, 'ArrayAggregator' );

    my $it = $self->select_iter(<<"SELECT");
select aggregate($on),$by from gff a
where exists (select count(*) from gff b where a.$by=b.$by)
and a.feature='$feature' group by $by
SELECT

    return sub {
        my $aggregate = $it->() or return;
        return ArrayAggregator->post_process($aggregate, $by, $on);
    };
}


=head2 select

 my $rows_aref = $gffstore->select("select seqname from gff");

Run raw select statement against db, return an arrayref of hashrefs

=cut
sub select{
    my $self = shift;
    my $stmt = shift;
    my $iter = $self->select_iter($stmt,@_);
    my @accum;
    while (my $row = $iter->()){
        push @accum,$row;
    }
    return \@accum;
}

=head2 select_row

 my $row_fref = $gffstore->select_row("select seqname from gff");

Run raw select statement against db, return a single row as hashref

=cut
sub select_row{
    my $self = shift;
    my $stmt = shift;

    my $dbh = $self->{dbh};
    my $sth = $dbh->prepare($stmt);
    $sth->execute(@_) or die "can't execute statement";
    
    my $row = $sth->fetchrow_hashref();
    $sth->finish;
    return $row;
}

=head2 select_row

 my $row_aref = $gffstore->select_col("select distinct seqname from gff");

Run raw select statement against db, return an aref of first column

=cut
sub select_col{
    my $self = shift;
    my $stmt = shift;

    my $dbh = $self->{dbh};
    my $sth = $dbh->prepare($stmt);
    $sth->execute(@_) or die "can't execute statement";
    
    my @accum;
    while (my $row = $sth->fetchrow_arrayref()){
        next if ! defined $row->[0];
        push @accum,$row->[0];
    }
    $sth->finish();
    return \@accum;
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
    say STDERR $select_stmt if $self->{debug};

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

=head2 overlappers

 overlappers($seqname, $feature, $start, $end);
 overlappers($seqname, 0,        $start, $end);

returns rows (arrayref of hashrefs) that overlaps with given region. 
Make sure you've indexed on ['seqname', 'feature', 'start', 'end'] 
or ['seqname', 'start', 'end'] if you're running many.

=cut
sub overlappers{
    my $self = shift;
    my ($seqname,$feature,$start,$end) = @_;
    die "need seqname, start, and end" unless ($seqname and $start and $end);
    my $dbh = $self->{dbh};
    state $sth_with_f = $dbh->prepare("select * from gff where seqname = ? and feature = ? and start <= ? and end >= ?");
    state $sth_without_f = $dbh->prepare("select * from gff where seqname = ? and start <= ? and end >= ?");

    if ($feature){
        $sth_with_f->execute($seqname,$feature,$end,$start);
        my $rv = $sth_with_f->fetchall_arrayref({});
        $sth_with_f->finish();
        return $rv;
    } else{
        $sth_without_f->execute($seqname,$end,$start);
        my $rv = $sth_without_f->fetchall_arrayref({});
        $sth_without_f->finish();
        return $rv;
    }
}

sub seqnames{
    my $self = shift;
    return $self->select_col("select distinct seqname from gff");
}

sub features{
    my $self = shift;
    return $self->select_col("select distinct feature from gff");
}

sub dump_gff{
    my $self = shift;
    my $dbh = $self->{dbh};
    my $select = $dbh->prepare("select * from gff");
    $select->execute();
    while (my $row = $select->fetchrow_arrayref()){
        say gff_to_string $row;
    }
}

sub DESTROY{
    my $self = shift;
    $self->{dbh}->disconnect;
    
    if ($self->{tmp}) {
        unlink0($self->{tmpfh}, $self->{dbname});
    }
}

1;

#my $RTREE_ERROR =  <<'END';
#==================================================================
#
#Searching for overlaps requires the R*Tree feature built into DBD::SQLite.  Please downloading the following file:
#
#http://search.cpan.org/CPAN/authors/id/A/AD/ADAMK/DBD-SQLite-1.31.tar.gz
#
#and rebuild with
#
#perl Makefile.PL -DSQLITE_ENABLE_RTREE=1 -DSQLITE_ENABLE_COLUMN_METADATA=1 -DSQLITE_ENABLE_FTS3=1 -DSQLITE_ENABLE_FTS3_PARENTHESIS=1
#make
#make test
#make install
#
#alternatively, you can put the following contents in ~/.cpan/prefs/sqlite-rtree.yml and do an 
#"install ADAMK/DBD-SQLite-1.31.tar.gz" at the cpan commandline.
#
#---
#comment: |
#  Compile RTREE feature
#
#match:
#  distribution: '^ADAMK/DBD-SQLite-\d'
#pl:
#  args: ["DEFINE=' -DSQLITE_ENABLE_RTREE=1 -DSQLITE_ENABLE_COLUMN_METADATA=1 -DSQLITE_ENABLE_FTS3=1 -DSQLITE_ENABLE_FTS3_PARENTHESIS=1 '"]
#depends:
#  configure_requires:
#    DBI: 1.58
#  requires:
#    DBI: 1.58
#
#==================================================================
#END

=head1 NAME
 
DZLab::Tools::GFFStore - DBI/DBD::SQLite backend for gff data
 
=head1 VERSION
 
This documentation refers to DZLab::Tools::GFFStore version 0.0.1
 
=head1 SYNOPSIS
 
    use DZLab::Tools::GFFStore;
  
=head1 DESCRIPTION
 
<description>

=head1 SUBROUTINES/METHODS 

<exports>

=over

=item 

=back
 
 
=head1 BUGS AND LIMITATIONS
 
Probably.

=cut
