#!/usr/bin/perl

use strict;
use warnings;
use JFR::Fasta;

MAIN: {
    my $fa = $ARGV[0] or die "usage: $0 FASTA\n";
    my $fp = JFR::Fasta->new($fa);
    while (my $rec = $fp->get_record()) {
        $rec->{'def'} =~ s/C_rubr_no_aliens_pep_Gene.(\d+).*/Crubr.$1/;
        $rec->{'def'} =~ s/^>.*_pep_(.*)/>$1/;
        $rec->{'def'} =~ s/\./|/;
        $rec->{'def'} =~ s/_//;
        $rec->{'def'} =~ s/ //g;
        print "$rec->{'def'}\n$rec->{'seq'}\n";
    }
}
