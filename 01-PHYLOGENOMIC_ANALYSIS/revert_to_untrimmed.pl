#!/usr/bin/perl

use strict;
use warnings;
use JFR::Fasta;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

our $VERSION = 0.01;

MAIN: {
    my $rh_opts = process_options();
    my $rh_s_ids = get_ids($rh_opts->{'subset_fasta'});
    my $rh_f_ids = get_ids($rh_opts->{'full_fasta'});
    check_that_all_subset_are_in_full($rh_s_ids,$rh_f_ids,$rh_opts);
    print_subset_of_full($rh_s_ids,$rh_opts->{'full_fasta'});
}

sub check_that_all_subset_are_in_full {
    my $rh_s_ids = shift;
    my $rh_f_ids = shift;
    my $rh_opts  = shift;
    my $absent_flag = 0;
    foreach my $id (keys %{$rh_s_ids}) {
        unless ($rh_f_ids->{$id}) {
            warn "$id is not present in $rh_opts->{'full_fasta'}\n";
            $absent_flag++;
         }
    }
    die "exiting\n" if ($absent_flag);
}

sub print_subset_of_full {
    my $rh_ids  = shift;
    my $full_fa = shift;
    my $fp      = JFR::Fasta->new($full_fa);
    while (my $rec = $fp->get_record()) {
        if ($rh_ids->{$rec->{'def'}}) {
            print "$rec->{'def'}\n$rec->{'seq'}\n";
        }
    }
}

sub get_ids {
    my $fa = shift;
    my %ids = ();
    my $fp = JFR::Fasta->new($fa);
    while (my $rec = $fp->get_record()) {
        $ids{$rec->{'def'}}++;
    }
    return \%ids;
}

sub usage  {
    print "revert_to_untrimmed.pl --full_fasta=FULL_FASTA --subset_fasta=SUBSET_FASTA [--version] [--help]\n";
    exit;
}

sub process_options {
    my $rh_opts = {};
    my $opt_results = Getopt::Long::GetOptions(
                              "version" => \$rh_opts->{'version'},
                                 "help" => \$rh_opts->{'help'},
                       "subset_fasta=s" => \$rh_opts->{'subset_fasta'},
                         "full_fasta=s" => \$rh_opts->{'full_fasta'},
                           );
    die "$VERSION\n" if ($rh_opts->{'version'});
    pod2usage({-exitval => 0, -verbose => 2}) if($rh_opts->{'help'});

    usage() unless $rh_opts->{'subset_fasta'};
    usage() unless $rh_opts->{'full_fasta'};
    return $rh_opts;
}

__END__

=head1 NAME

B<revert_to_untrimmed.pl> - Revert to untrimmed

=head1 AUTHOR

Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

revert_to_untrimmed.pl --full_fasta=FULL_FASTA --subset_fasta=SUBSET_FASTA [--version] [--help]

=head1 DESCRIPTION

This program takes 2 fasta files. It grabs the ids from the SUBSET_FASTA file and then prints the corresponding sequences from the FULL_FASTA. All of the ids in the SUBSET_FASTA file should be in the FULL_FASTA file.

=head1 BUGS

Please report them to the author.

=head1 COPYRIGHT

Copyright (C) 2020, Joseph F. Ryan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut
