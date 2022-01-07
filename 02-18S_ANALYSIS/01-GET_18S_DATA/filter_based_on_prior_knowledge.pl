#!/usr/bin/perl

# script removes several problematic 18S Sequences.
# justification is provided in the manuscript

use strict;
use warnings;
use JFR::Fasta;

# manually adds Nematostella vectensis (XR_004291954.1) at the end

MAIN: {
    my $file = $ARGV[0] or die "usage: $0 FASTA\n";
    my $fp = JFR::Fasta->new($file);
    while (my $rec = $fp->get_record()) {
        next if ($rec->{'def'} eq '>Virgularia_gustaviana_OC');
        next if ($rec->{'def'} eq '>Alatina_philippina_CU');
        next if ($rec->{'def'} eq '>Junceella_aquamata_OC');
        next if ($rec->{'def'} eq '>Junceella_fragilis_OC');
        next if ($rec->{'def'} eq '>Subergorgia_ornata_OC');
        if ($rec->{'def'} eq '>Carybdea_marsupialis_CU') {
            print ">Alatinidae_indet._CU\n$rec->{'seq'}\n";
        } elsif ($rec->{'def'} eq '>Darwin_sp_CU') {
            print ">Gerongia_rifkinae_CU\n$rec->{'seq'}\n";
        } else {
            print "$rec->{'def'}\n$rec->{'seq'}\n";
        }
    }
    print ">Nematostella_vectensis_HE\ntatctggttgatcctgccagtagtcatatgcttgtctcaaagattaagccatgcatgtctaagtataagcacttgtactgtgaaactgcgaatggctcattaaatcagttatcgtttatttgattgtacctttactacttggataaccgtggtaattctagagctaatacatgcgaaaagtcccgacttctggaagggatgtatttattagattcaaaaccaatgcgggttctgcccggtcatttggtgattcatagtaactgttcgaatcgcagggcctcgcgccggcgatgtttcattcaaatttctgccctatcaactgtcgatggtaaggtgttggcttaccatggttacaacgggtgacggagaattagggttcgattccggagagggagcctgcgaaacggctatcacatccaaggaaggcagcaggcgcgcaaattacccaatcctgactcagggaggtagtgacaagaaataacaatacagggcttttctaagtcttgtaattggaatgagtacaacttaaatcctttaacgaggatccattggagggcaagtctggtgccagcagccgcggtaattccacctccaatagcgtatattaaagttgttgcagttaaaaagctcgtagttggatttcgggacggcacggtcggtccaccgcaaggtgtgtcactggccgggctgttcttcctcgcaaagactgcgtgtgctcttagctgagtgtgcgcaggacttgcgacgtttactttgaaaaaattagagtgttcaaagcaggccatcgcttgaatacataagcatggaataatggaataggactttggttctattttgttggtttctggaaccgaagtaatgattaagagggacagttgggggcattcgtatttcgttgtcagaggtgaaattcttggatttacgaaagacgaactactgcgaaagcatttgccaagaatgttttcattaatcaagaacgaaagttagaggatcgaagacgatcagataccgtcctagttctaaccataaacgatgccgactagggatcagagagtgttattggatgacctctttggcaccttatgggaaaccaaagtttttgggttccgggggaagtatggttgcaaagctgaaacttaaaggaattgacggaagggcaccaccaggagtggagcctgcggcttaatttgactcaacacggggaaactcaccaggtccagacataggaaggattgacagattgagagctctttcttgattctatgggtggtggtgcatggccgttcttagttggtggagtgatttgtctggttaattccgttaacgaacgagaccttaacctgctaaatagttacgctaatctcgattggcggctaacttcttagagggactgttggtgttcaaccaaagtcaggaaggcaataacaggtctgtgatgcccttagatgttctgggccgcacgcgcgctacactgacggtgtcaacgagtctttccttcgccggaaggcgtgggtaatcttgtgaaatatcgtcgtgctggggatagatcattgcaattattgatcttgaacgaggaattcctagtaagcgcgagtcatcagctcgcgttgattacgtccctgccctttgtacacaccgcccgtcgctactaccgattgaatggtttagtgaggccttctgattggcgccgcggccccggcaacggagccacggattgtcgaaaagttggtcaaacttgatcatttagaggaagtaaaagtcgtaacaaggtttccgtaggtgaacctgcggaaggatcatta\n";


}
