#!/usr/bin/env perl
# ___UNDOCUMENTED___

use warnings;
use strict;
use Data::Dumper;
use Carp;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/sum/;

my $DATA_HANDLE = 'ARGV';
my $output;
my $append; # default behaviour is to merge
my $dnf;

# Grabs and parses command line options
my $result = GetOptions (
    'dnf|d'       => \$dnf,
    'output|o=s'  => \$output,
    'verbose|v'   => sub { use diagnostics; },
    'quiet|q'     => sub { no warnings; },
    'help|h'      => sub { pod2usage ( -verbose => 1 ); },
    'manual|m'    => sub { pod2usage ( -verbose => 2 ); }
);

# Check required command line parameters
pod2usage ( -verbose => 1 )
unless @ARGV and $result;

if ($output) {
    open my $USER_OUT, '>', $output or croak "Can't open $output for writing: $!";
    select $USER_OUT;
}

my %data;
my @headers;
my @fields;

while (<$DATA_HANDLE>) {    
    unless (m/\d+/) {
        @headers = split /\t/;
        next;
    }
    else {
        @fields = split /\t/;
    }

    croak "Number of fields: ", scalar @fields, " is different from number of headers: ", scalar @headers, "\n",
    join "\t", @headers, "\n", join "\t", @fields
    if scalar @headers and @fields and  @fields != scalar @headers;

    for my $stat_index (0 .. @headers - 1) {
        push @{$data{$headers[$stat_index]}}, $fields[$stat_index];
    }

}
my $processed_data = $dnf ? process_dnf_fields (\%data) : process_tnf_fields (\%data);

my @sorted_keys = sort keys %{$processed_data};
my $headers = join "\t", @sorted_keys;
my $fields  = join "\t", map { $processed_data->{$_} } @sorted_keys;

print $headers, "\n", $fields, "\n";



sub process_tnf_fields {
    my $data = shift;
    my %freqs;

    $freqs{bp}     = sum @{$data{'bp'}};

    $freqs{cg}     = sum @{$data{'CG'}};
    $freqs{chg}    = sum @{$data{'CHG'}};
    $freqs{chh}    = sum @{$data{'CHH'}};
    $freqs{c}      = sum ($freqs{cg}, $freqs{chg}, $freqs{chh});

    $freqs{tg}     = sum @{$data{'TG'}};
    $freqs{thg}    = sum @{$data{'THG'}};
    $freqs{thh}    = sum @{$data{'THH'}};
    $freqs{t}      = sum ($freqs{tg}, $freqs{thg}, $freqs{thh});

    $freqs{f_cg}    = sum @{$data{'filtered_CG'}};
    $freqs{f_chg}    = sum @{$data{'filtered_CHG'}};
    $freqs{f_chh}    = sum @{$data{'filtered_CHH'}};
    $freqs{f_c}     = sum ($freqs{f_cg}, $freqs{f_chg}, $freqs{f_chh});

    $freqs{f_tg}    = sum @{$data{'filtered_TG'}};
    $freqs{f_thg}    = sum @{$data{'filtered_THG'}};
    $freqs{f_thh}    = sum @{$data{'filtered_THH'}};
    $freqs{f_t}     = sum ($freqs{f_tg}, $freqs{f_thg}, $freqs{f_thh});

    $freqs{c_ratio}  = $freqs{c} / ($freqs{c} + $freqs{t});
    $freqs{cg_ratio} = $freqs{cg} / ($freqs{cg} + $freqs{tg});
    $freqs{chg_ratio} = $freqs{chg} / ($freqs{chg} + $freqs{thg});
    $freqs{chh_ratio} = $freqs{chh} / ($freqs{chh} + $freqs{thh});

    $freqs{f_c_ratio}  = $freqs{f_c} / ($freqs{f_c} + $freqs{f_t});
    $freqs{f_cg_ratio} = $freqs{f_cg} / ($freqs{f_cg} + $freqs{f_tg});
    $freqs{f_chg_ratio} = $freqs{f_chg} / ($freqs{f_chg} + $freqs{f_thg});
    $freqs{f_chh_ratio} = $freqs{f_chh} / ($freqs{f_chh} + $freqs{f_thh});

    return \%freqs;
}



sub process_dnf_fields {
    my $data = shift;
    my %freqs;

    $freqs{bp}    = sum @{$data{'bp'}};

    $freqs{cg}    = sum @{$data{'CG'}};
    $freqs{ca}    = sum @{$data{'CA'}};
    $freqs{cc}    = sum @{$data{'CC'}};
    $freqs{ct}    = sum @{$data{'CT'}};
    $freqs{c}     = sum ($freqs{ca}, $freqs{cc}, $freqs{cg}, $freqs{ct});
    $freqs{no_cg} = sum ($freqs{ca}, $freqs{cc}, $freqs{ct});

    $freqs{tg}    = sum @{$data{'TG'}};
    $freqs{ta}    = sum @{$data{'TA'}};
    $freqs{tc}    = sum @{$data{'TC'}};
    $freqs{tt}    = sum @{$data{'TT'}};
    $freqs{t}     = sum ($freqs{ta}, $freqs{tc}, $freqs{tg}, $freqs{tt});
    $freqs{no_tg} = sum ($freqs{ta}, $freqs{tc}, $freqs{tt});

    $freqs{f_cg}    = sum @{$data{'filtered_CG'}};
    $freqs{f_ca}    = sum @{$data{'filtered_CA'}};
    $freqs{f_cc}    = sum @{$data{'filtered_CC'}};
    $freqs{f_ct}    = sum @{$data{'filtered_CT'}};
    $freqs{f_c}     = sum ($freqs{f_ca}, $freqs{f_cc}, $freqs{f_cg}, $freqs{f_ct});
    $freqs{f_no_cg} = sum ($freqs{f_ca}, $freqs{f_cc}, $freqs{f_ct});

    $freqs{f_tg}    = sum @{$data{'filtered_TG'}};
    $freqs{f_ta}    = sum @{$data{'filtered_TA'}};
    $freqs{f_tc}    = sum @{$data{'filtered_TC'}};
    $freqs{f_tt}    = sum @{$data{'filtered_TT'}};
    $freqs{f_t}     = sum ($freqs{f_ta}, $freqs{f_tc}, $freqs{f_tg}, $freqs{f_tt});
    $freqs{f_no_tg} = sum ($freqs{f_ta}, $freqs{f_tc}, $freqs{f_tt});

    $freqs{c_ratio}  = $freqs{c} / ($freqs{c} + $freqs{t});
    $freqs{ca_ratio} = $freqs{ca} / ($freqs{ca} + $freqs{ta});
    $freqs{cc_ratio} = $freqs{cc} / ($freqs{cc} + $freqs{tc});
    $freqs{cg_ratio} = $freqs{cg} / ($freqs{cg} + $freqs{tg});
    $freqs{ct_ratio} = $freqs{ct} / ($freqs{ct} + $freqs{tt});

    $freqs{f_c_ratio}  = $freqs{f_c} / ($freqs{f_c} + $freqs{f_t});
    $freqs{f_ca_ratio} = $freqs{f_ca} / ($freqs{f_ca} + $freqs{f_ta});
    $freqs{f_cc_ratio} = $freqs{f_cc} / ($freqs{f_cc} + $freqs{f_tc});
    $freqs{f_cg_ratio} = $freqs{f_cg} / ($freqs{f_cg} + $freqs{f_tg});
    $freqs{f_ct_ratio} = $freqs{f_ct} / ($freqs{f_ct} + $freqs{f_tt});

    return \%freqs;
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

 $Rev: $:
 $Author: $:
 $Date: $:
 $HeadURL: $:
 $Id: $:

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
