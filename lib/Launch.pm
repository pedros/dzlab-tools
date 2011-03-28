package Launch;
use version; our $VERSION = qv('0.0.1');
use strict;
use warnings;
use Data::Dumper;
use feature 'say';
use Carp;
use Cwd;
use File::Spec::Functions;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use IPC::Open3;
use Log::Log4perl qw/get_logger/;
use Parallel::ForkManager;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw();
our @EXPORT = qw(launch plaunch);

=head2 launch

 expected - This is a file (or arrayref of files) which we expect to be produce.
 force    - Run even if there is valid hash in donedir and file exists.
 donedir  - The place to put. default to cwd.

=cut

sub launch{
    my ($cmd, %opt) = @_;
    my $logger    = get_logger("Launch");

    my $hash      = md5_hex($cmd);
    my $done_dir  = delete $opt{donedir} // 0;
    my $done_file = $done_dir ? catfile($done_dir,$hash) : '';
    my $force     = delete $opt{force} // 0;

    $logger->info("running [$cmd], done = $done_file, force = $force");

    my @expected;
    if (exists $opt{expected}){
        if (ref $opt{expected} eq 'ARRAY'){ 
            @expected = @{$opt{expected}}
        } else {
            @expected = ($opt{expected});
        }
        delete $opt{expected}
    }

    die "unknown parameters passed to doit" . Dumper \%opt if (%opt);

    if (!$force){
        if ((! $done_dir || -f $done_file) && (! scalar(@expected) || grep {-f} @expected)){
            $logger->info("Already done, skipping: '$cmd' ");
            return 1;
        }
    }
    
    if (0==system($cmd)){
        mkdir $done_dir if ($done_dir && ! -d $done_dir); 
        
        if (! @expected || grep {-f} @expected){ 
            if ($done_dir){
                open my $out, '>', $done_file 
                    or $logger->logdie("finished, but can't open $done_file");
                print $out $cmd;
                close $out
                    or $logger->logdie("finished, but can't close $done_file");
            }
            $logger->info("Successfully launched and finished [$cmd]");
        } else {
            $logger->logdie("command seems to have run but expected files not produced [$cmd]");
        }
    } else {
        $logger->logdie("failed to run, dying: ($cmd]");
    }
}
sub plaunch{
    my ($numprocs, @jobs) = @_;
    # Max 30 processes for parallel download
    my $pm = new Parallel::ForkManager($numprocs);

    foreach my $j (@jobs) {
        $pm->start and next; # do the fork
        launch(@$j);
        $pm->finish; # do the exit in the child process
    }
    $pm->wait_all_children;
}
1;
