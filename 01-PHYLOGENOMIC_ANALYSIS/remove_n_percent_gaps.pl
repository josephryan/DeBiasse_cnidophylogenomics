#!/usr/bin/perl

use strict;
use warnings;
use JFR::Fasta;

our $GAP = '-';

MAIN: {
    my $fa = $ARGV[0] or die "usage: $0 FASTA PERCENT_GAPS_ALLOWED\n";
    my $per = $ARGV[1] or die "usage: $0 FASTA PERCENT_GAPS_ALLOWED\n";
    remove_n_percent_gaps($fa,$per);
}

sub remove_n_percent_gaps {
    my $fa = shift;
    my $per = shift;
    $per /= 100;
    my $fp = JFR::Fasta->new($fa);
    while (my $rec = $fp->get_record()) {
        my $count = () = $rec->{'seq'} =~ /\Q$GAP/g;
        my $percent_gaps = $count / length($rec->{'seq'});
        print "$rec->{'def'}\n$rec->{'seq'}\n" if ($percent_gaps <= $per);
    }
}
