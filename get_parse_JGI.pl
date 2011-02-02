#!/usr/bin/env perl
# ___UNDOCUMENTED___

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use LWP::Simple qw(get);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV;

my $fid;
my $from = 1;
my $db;
my $table;
my $output;

# Grabs and parses command line options
my $result = GetOptions (
    'fid|i=s'    => \$fid,
    'db|d=s'     => \$db,
    'from|f=i'   => \$from,
    'table|t=s'  => \$table,
    'output|o=s' => \$output,
    'verbose|v'  => sub { use diagnostics; },
    'quiet|q'    => sub { no warnings; },
    'help|h'     => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'   => sub { pod2usage ( -verbose => 2 ); }
);

if ($output) {
    open my $USER_OUT, '>', $output or carp "Can't open $output for writing: $!";
    select $USER_OUT;
}

if ($fid) {

    if (-e $fid) {
        open my $FID, '<', $fid or croak "Can't open $fid: $!";
        while (<$FID>) {
            chomp;
            my $target
            = 'http://genome.jgi-psf.org/cgi-bin/colorSeqViewer?db=' . $db . '&table=' . $table . '&fid=' . $_;
            get_parse ($target, $_);
        }
    }
    else {
        my $target
        = 'http://genome.jgi-psf.org/cgi-bin/colorSeqViewer?db=' . $db . '&table=' . $table . '&fid=' . $fid;
        get_parse ($target, $fid);
    }
}
elsif (@ARGV) {
    for ($ARGV[0] .. $ARGV[$#ARGV]) {
        my $target
        = 'http://genome.jgi-psf.org/cgi-bin/colorSeqViewer?db=' . $db . '&table=' . $table . '&fid=' . $_;
        get_parse ($target, $_);
    }
}
else {
    $fid = $from;
    while (1) {
        print STDERR "Fetching $fid\n" unless $fid % 100;
        my $target
        = 'http://genome.jgi-psf.org/cgi-bin/colorSeqViewer?db=' . $db . '&table=' . $table . '&fid=' . $fid;
        get_parse ($target, $fid++);
    }
}


sub get_parse {
    my ($target, $fid, $html) = @_;

    $html = get ( $target ) until $html;

    return 0 if $html =~ m/Software error/;
    
    my @targets = grep { /^name|^strand|^scaff/ } (split /[\n\r]/, $html);

    print join "\t", $fid, 'JGI', 'repeat', '';

    for (@targets[2, 1, 3, 3, 0]) {
        chomp;
        s/<br>//;
        
        if (/\(padded\)/) {
            print q{.}, "\t";
            next;
        }
        
        my ($id, $value) = split /:/;

        print "name=" if $id =~ /name/;

        if ($id =~ m/scaffStart\/end/) {
            my ($start, $end) = split /,/, $value;
            print join "\t", $start, $end, '';
        }
        else {
            print $value, "\t";
        }
    }
    print "\n";
    return 1;
}


__END__


=head1 NAME

 name.pl - Short description

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 OPTIONS

 name.pl [OPTION]... [FILE]...

 -o, --output      filename to write results to (defaults to STDOUT)
 -v, --verbose     output perl's diagnostic and warning messages
 -q, --quiet       supress perl's diagnostic and warning messages
 -h, --help        print this information
 -m, --manual      print the plain old documentation page

=head1 REVISION

 Version 0.0.1

 $Rev: 249 $:
 $Author: psilva $:
 $Date: 2010-01-11 21:24:34 -0800 (Mon, 11 Jan 2010) $:
 $HeadURL: http://dzlab.pmb.berkeley.edu/svn/bisulfite/trunk/get_parse_JGI.pl $:
 $Id: get_parse_JGI.pl 249 2010-01-12 05:24:34Z psilva $:

=head1 AUTHOR

 Pedro Silva <psilva@nature.berkeley.edu/>
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
