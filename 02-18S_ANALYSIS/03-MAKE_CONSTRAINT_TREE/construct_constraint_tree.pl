#!/usr/bin/perl

use strict;
use warnings;
use JFR::Fasta;
use Data::Dumper;

our $VERSION = 0.03;

our $TRANSCRIPT_CSV = 'cnidophylo_names_abbrev.csv';
our $SSRNA_FA = '../02-ALN/ssu.dir.eukarya.mask.stk.fa';
our $TRANSCRIPT_TRE = '../../01-PHYLOGENOMIC_ANALYSIS/02-TREES/cnido_748_untrimmed_b.tre';

our %OUTGROUPS = ('Ctele' => 1, 'Lgiga' => 1, 'Dmela' => 1, 'Spurp' => 1,
                  'Tgutt' => 1, 'Tadha' => 1, 'Mleid' => 1, 'Aquee' => 1);

our %MANUAL = ('Cassiopeasp' => 'Cassiopea_frondosa_SC',
               'Turritopsissp' => 'Turritopsis_dohrnii_HY',
               'Renilla' => 'Renilla_reniformis_OC',
               'Bcave' => 'Bunodosoma_grande_HE',
               'Bmexi' => 'Botruanthus_benedeni_CE',
               'Caust' => 'Corynactis_californica_HE',
               'Eflex' => 'Eunicea_laciniata_OC',
               'Fliza' => 'Favia_fragum_HE',
               'Hdigi' => 'Hormathia_armata_HE',
               'Leptog' => 'Leptogorgia_chilensis_OC',
               'Lnegl' => 'Lebrunia_coralligens_HE',
               'Paste' => 'Porites_cylindrica_HE',
               'Plcarn' => 'Platygyra_daedalea_HE',
               'Rindo' => 'Rhodactis_rhodostoma_HE',
               'Ryuma' => 'Ricordea_florida_HE',
               'Cbore' => 'Pachycerianthus_borealis_CE',
               );


#                  'Hpoly' => 'Hydractinia polyclina',
our %AMBIGUOUS = (
                  'Adigi' => 'Acropora_abrolhosensis_HE',
                  'Caust' => 'Corynactis_australis_HE',
                  'Cbras' => 'Ceriantheomorphe_brasiliensis_CE',
                  'Ccapi' => 'Cyanea_capillata_SC',
                  'Chemi' => 'Clytia_hemisphaerica_HY',
                  'Cpoly' => 'Calliactis_polypus_HE',
                  'Equad' => 'Entacmaea_quadricolor_HE',
                  'Maure' => 'Madracis_kirbyi_HE',
                  'Pvari' => 'Palythoa_variabilis_HE',
                  'Phydr' => 'Polypodium_hydriforme_PO',
                  'Seleg' => 'Sagartia_elegans_HE',
                  );

# Hpoly represents Hydractinia polyclina which is not in our 18s set
#     but Hpoly erroneously matches Hydra_polymorphus_HY
# Crubr represents Corallium rubrum which is not in our 18S set
#     (existing sequences are 626 and 627bp
our %SKIP = ('Hpoly' => 1, 'Crubr' => 1);

# These are skipped because they were removed from the phylogenomics
#     analysis because of low BUSCO scores or contamination
our %SKIP2 = ('Pperi' => 1,
             'Hcris' => 1,
             'Scocc' => 1,
             'Aviri' => 1,
             'Mfave' => 1,
             'Muri'  => 1,
             'Hcoer' => 1,
            );

MAIN: {
    # { 'Tkita' => 'Thelohanellus kitauei', ...}
    my $rh_tp = get_transcript_taxa($TRANSCRIPT_CSV);

    # { 'Pseudocrypthelia_pachypoma_HY' => 1, ...}
    my $rh_ss = get_ssrna_taxa($SSRNA_FA);

    # { 'Hvulg' => 'Hydra_vulgaris_HY', ...}
    my $rh_match = get_matches($rh_tp,$rh_ss);

    add_manual($rh_match);

    prune_missing_and_outgroup($rh_tp,$rh_match,$TRANSCRIPT_TRE);
    my $tre_str = remove_support_and_branch_lengths("$TRANSCRIPT_TRE.pruned");
    my $constraint = change_names($tre_str,$rh_match);
    print $constraint . ';' . "\n";
}

sub change_names {
    my $tre  = shift;
    my $rh_m = shift;
    my %rev  = ();  
    my $tmp_tre = $tre;
    $tmp_tre =~ s/[\(\)\,]+/,/g;
    my @taxa = split /,/, $tmp_tre;
    foreach my $taxon (@taxa) {
        next if ($taxon =~ m/^\s*$/);
        my $num = get_num_matches($tre,$taxon);
        die "unexpected match# $taxon:$num" unless ($num == 1);        
        die "unexpected absence from \$rh_m:$taxon" unless ($rh_m->{$taxon});
        $tre =~ s/$taxon/$rh_m->{$taxon}/;
    }
    return $tre;
}

sub get_num_matches {
    my $str = shift;
    my $regex = shift;
    my $num_matches = 0;
    while ($str =~ m/$regex/g) {
        $num_matches++;
    }
    return $num_matches;
}

sub add_manual {
    my $rh_match = shift;
    foreach my $key (keys %MANUAL) {
        if ($rh_match->{$key}) {
            die "unexpected $key: $MANUAL{$key} $rh_match->{$key}" if ($rh_match->{$key});
        }
        $rh_match->{$key} = $MANUAL{$key};
    }
}

sub remove_support_and_branch_lengths {
    my $tre = shift;
    my $str = '';
    open IN, $tre or die "cannot open $tre:$!";
    while (my $line = <IN>) {
        chomp $line;
        $str .= $line;
    }
    unlink $tre;
    # these regexes are borrowed from AfterPhylo version 0.9.1 
    #    https://github.com/qiyunzhu/AfterPhylo
    $str =~ s/\[.+?\]//g;
    $str =~ s/\)[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?([:,\)\[])/\)$2/g;
    $str =~ s/:[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?([,)[])/$2/g;
    $str =~ s/;//g; # not from AfterPhylo
    return $str;
}

sub prune_missing_and_outgroup {
    my $rh_tp = shift;
    my $rh_m  = shift;
    my $tree  = shift;
    my $count = 666;
    my @to_remove = ();
    foreach my $key (sort keys %{$rh_tp}) {
        next if ($rh_m->{$key});
        next if ($SKIP2{$key});
        push @to_remove, $key;
    }
    my $string = join ' ', @to_remove;
    my $cmd = "/usr/local/phyutility/phyutility -pr -in $tree -out $tree.pruned -names $string";
    system $cmd;
}

sub get_matches {
    my $rh_tp = shift;
    my $rh_ss = shift;
    my %matches = ();

    my %tp_vals = ();
    foreach my $key (keys %{$rh_tp}) {
        $rh_tp->{$key} =~ s/ /_/g;
        $tp_vals{$rh_tp->{$key}} = $key;
    }

    foreach my $ss (sort keys %{$rh_ss}) {
        die "unexpected mismatch:$ss" unless $ss =~ m/^(.*)_[A-Z][A-Z]$/;
        my $name  = $1;
        # done separately because not all have 4-letter sp. 
        die "unexpected" unless ($ss =~ m/^([^_])[^_]+_(.{4})/);
        my $abrev = $1 . $2;

        next if ($SKIP{$abrev});
        next if ($SKIP2{$abrev});
        if ($MANUAL{$abrev}) {
            $matches{$abrev} = $MANUAL{$abrev};
        } elsif ($tp_vals{$name}) { 
            $matches{$tp_vals{$name}} = $ss;
        } elsif ($AMBIGUOUS{$abrev}) {
            $matches{$abrev} = $AMBIGUOUS{$abrev};
        } elsif ($rh_tp->{$abrev}) {
            $matches{$abrev} = $ss;
        }
    }
    return \%matches;
}

sub get_ssrna_taxa {
    my $fa = shift;
    my %ss = ();
    my $fp = JFR::Fasta->new($fa);
    while (my $rec = $fp->get_record()) {
        my $id = JFR::Fasta->get_def_w_o_gt($rec->{'def'});
        $ss{$id}++;
    }
    return \%ss;
}

sub get_transcript_taxa {
    my $tfa = shift;
    my %taxa = ();
    open IN, $tfa or die "cannot open $tfa:$!";
    while (my $line = <IN>) {
        chomp $line; 
        my @fields = split /,/, $line;
        $taxa{$fields[0]} = $fields[1];
    }
    return \%taxa;
}
