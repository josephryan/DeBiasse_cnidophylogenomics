#!/usr/bin/perl

# 18s.long.fa will contain more information about sequence

use strict;
use FileHandle;
use IO::Uncompress::Gunzip qw($GunzipError);
use warnings;
use Data::Dumper;

our $VERSION = 0.1;

our $FILE = '18s.cnidaria.genbank.gz';
our $REQUIRE_CNID_CLASS = 1;
our %SKIP_ACC = ('AY935208' => 1, 'AF099104' => 1, 'MF322725' => 1);
our $MAX_AMBIG  = 5;
our $MIN_LEN    = 1000;
our $OUTFILE = '18s.long.fa';
our $REJECTS = 'rejected_because_of_deflines.txt';

MAIN: {
    my $rh_data = get_data($FILE);
    my $rh_final = finalize_data($rh_data);    
    print_fasta($rh_final);
}

sub get_genera {
    my $rh_d   = shift;
    my %genera = ();
    foreach my $key (keys %{$rh_d}) {
        my @fields = split /\s+/, $key;
        $genera{$fields[0]}++ if (scalar(@fields) == 2);   
    }
    return \%genera;
}

sub finalize_data {
    my $rh_d  = shift;
    my %final = ();
    foreach my $org (keys %{$rh_d}) {
        foreach my $acc (keys %{$rh_d->{$org}}) {
            my $start = 1;
            my $len   = $rh_d->{$org}->{$acc}->{'len'};
            my $end   = $len;
            my $seq   = $rh_d->{$org}->{$acc}->{'seq'};
#            my $end   = length($rh_d->{$org}->{$acc}->{'seq'});
            foreach my $rh_r (@{$rh_d->{$org}->{$acc}->{'rrna'}}) {
                if ($rh_r->{'misc'} =~ m/18[Ss]/) {
                    $start = $rh_r->{'start'};
                    $end   = $rh_r->{'end'};
                    $len   = $rh_r->{'end'} - $start + 1;
                    $seq   = substr($rh_d->{$org}->{$acc}->{'seq'},($start - 1),$len);
                }
            }
#ZZZ
#            $end   = $end   || length($rh_d->{$org}->{$acc}->{'seq'});
#            my $len = $end - $start;
#            my $seq   = substr($rh_d->{$org}->{$acc}->{'seq'},($start - 1),$len);
            $final{$org}->{$acc}->{'seq'} = $seq;
            $final{$org}->{$acc}->{'len'} = $len;
            $final{$org}->{$acc}->{'tax'} = $rh_d->{$org}->{$acc}->{'tax'};
        }
    }
    return \%final;
}

sub print_fasta {
    my $rh_d = shift;
    my %printed = ();
    my $rh_g = get_genera($rh_d);
    open OUT, ">$OUTFILE" or die "cannot open >$OUTFILE:$!";
    foreach my $org (sort keys %{$rh_d}) {
        next if ($org =~ m/(environ)|(parasite)|(proliferative)/);
        foreach my $acc (sort {$rh_d->{$org}->{$b}->{'len'} <=>
                               $rh_d->{$org}->{$a} } keys %{$rh_d->{$org}}) {
            next if ($printed{$org});
            my $len = $rh_d->{$org}->{$acc}->{'len'};
#ZZZ
my $start = $rh_d->{$org}->{$acc}->{'start'};
my $end = $rh_d->{$org}->{$acc}->{'end'};
die "$org $acc $len $start $end" unless $len;
            next if ($len < $MIN_LEN);
            my $tax = get_tax($rh_d,$org,$acc);
            next if ($REQUIRE_CNID_CLASS && !$tax);
            my @fields = split /\s+/, $org;
            my $short_org = "$fields[0]_$fields[1]";
            if (scalar(@fields) == 2) {
                if ($fields[1] eq 'sp.' && $rh_g->{$fields[0]}) {
                    next;
                } elsif ($fields[1] eq 'sp.' && !$rh_g->{$fields[0]}) {
                    $fields[1] =~ s/\.$//;
                    $short_org = "$fields[0]_$fields[1]"; 
                    $rh_g->{$fields[0]}++;
                }
            } else {
                if ($fields[1] eq 'sp.' && $rh_g->{$fields[0]}) {
                    next;
                } elsif ($fields[1] eq 'sp.' && !$rh_g->{$fields[0]}) {
                    $fields[1] =~ s/\.$//;
                    $short_org = "$fields[0]_$fields[1]"; 
                    $rh_g->{$fields[0]}++;
                } elsif ($fields[1] eq 'cf.' && $rh_g->{$fields[0]}) {
                    next;
                } elsif ($fields[1] eq 'cf.' && !$rh_g->{$fields[0]}) {
                    $fields[1] =~ s/\.$//;
                    $short_org = "$fields[0]_$fields[1]"; 
                    $rh_g->{$fields[0]}++;
                } elsif ($fields[0] eq 'cf.') {
                    if ($rh_g->{$fields[1]}) {
                        next;
                    } else {
                        $fields[0] =~ s/\.$//;
                        $short_org = "$fields[1]_$fields[0]";
                        $rh_g->{$fields[1]}++;
                    }
                } else {
                    if ($rh_g->{$fields[0]}) {
                        next;
                    } else {
                        $short_org = "$fields[0]_sp"; 
                        $rh_g->{$fields[0]}++;
                    }
                }
            }


            my $name = $org . ".$acc";
            $name =~ s/'//g;
            $short_org =~ s/'//g;
            $name =~ s/ /_/g;
            my $def = "$name ($tax)";
            my $class = get_class_code($tax);
            print OUT ">$def\n$rh_d->{$org}->{$acc}->{'seq'}\n";
            print ">${short_org}_$class\n$rh_d->{$org}->{$acc}->{'seq'}\n";
            $printed{$org}++;
        }
    }
}

sub get_class_code {
    my $tax = shift;
    if ($tax =~ m/staurozoa/i) {
        return "ST";
    } elsif ($tax =~ m/cubozoa/i) {
        return "CU";
    } elsif ($tax =~ m/scyphozoa/i) {
        return "SC";
    } elsif ($tax =~ m/hydrozoa/i) {
        return "HY";
    } elsif ($tax =~ m/hexacorallia/i) {
        if ($tax =~ m/ceriantharia/i) {
            return "CE";
        } else {
            return "HE";
        }
    } elsif ($tax =~ m/octocorallia/i) {
        return "OC";
    } elsif ($tax =~ m/myxozoa/i) {
        return "MY";
    } elsif ($tax =~ m/polypodiozoa/i) {
        return "PO";
    } elsif ($tax =~ m/unspecified/i) { # Relicanthus daphneae
        return "UN";
    } else {
        die "unexpected $tax";
    }
} 

sub get_tax {
    my $rh_d = shift;
    my $org = shift;
    my $acc = shift;

    my $tax = '';
    if ($org eq 'Relicanthus daphneae') {
        $tax = 'unspecified';
    } elsif (!$rh_d->{$org}->{$acc}->{'tax'}->[3]) {
        return '';
    } else {
        $tax = "$rh_d->{$org}->{$acc}->{'tax'}->[3]";
        $tax .= "/$rh_d->{$org}->{$acc}->{'tax'}->[4]" if ($rh_d->{$org}->{$acc}->{'tax'}->[4]);
        $tax .= "/$rh_d->{$org}->{$acc}->{'tax'}->[5]" if ($rh_d->{$org}->{$acc}->{'tax'}->[5]);
    }
    return $tax;
}

sub get_data {
    my $file = shift;
    my $fh = IO::Uncompress::Gunzip->new($file) or die "IO::Uncompress::Gunzip of $file failed: $GunzipError\n";
    my %data = ();
    my @cur = ();
    open(my $fho, ">", "$REJECTS") or die "cannot open >$REJECTS:#!";
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ m/^LOCUS/) {
            process_cur(\@cur,\%data,$fho) if scalar(@cur);
            @cur = ();
        }
        push @cur, $line;
    }
    process_cur(\@cur,\%data,$fho);
    return \%data;
}
sub process_cur {
    my $ra_lines = shift;
    my $rh_data  = shift;
    my $fh = shift;
    my $first = shift @{$ra_lines};
    my @f = split /\s+/, $first;
    my $acc = $f[1];
    return if ($SKIP_ACC{$acc});
    my $len = $f[2];
#    die "unexpected: $first" unless ($f[4] eq 'DNA');
#    die "unexpected: $first" unless ($f[5] eq 'linear');
    my $source = '';
    my $org    = '';
    my $seq_flag = 0;
    my $seq = '';
    my $def_flag = 0;
    my $def = '';
    my $class = '';
    my $tax_flag = 0;
    my $tax_line = '';
    my @tax = ();
    my $hit_ref = 0;
    my $feat_flag = 0;
    my %features;
    my @rrna = ();
    my $rrna_count = -1;
    foreach my $line (@{$ra_lines}) {
        if ($line =~ m/^\s+ORGANISM\s+(\S+.*)/) {
            next if ($org);
            return 0 if ($org =~ m/environmental/);
            $org      = $1;
            $hit_ref  = 0;
            $tax_flag = 1;
        } elsif ($line =~ m/^ACCESSION/) {
            if ($def =~ m/(28S)|(spacer)|(NADH)/i) {
                print $fh "28:$def\n";
                $def = '';
                return 0;
            } elsif ($def =~ m/(18S)|(small subunit ribosomal)/i) {
                $def = '';
            } else {
                print $fh "$def\n";
                $def = '';
                return 0;
            }
        } elsif ($line =~ m/^DEFINITION/) {
            $def .= $line;
        } elsif ($def) {
            $line =~ s/^\s+/ /;
            $def .= $line;
        } elsif ($line =~ m/^ORIGIN/) {
            $seq_flag  = 1;
            $feat_flag = 0;
        } elsif ($line =~ m/^CONTIG/) {
            $feat_flag = 0;
        } elsif ($line =~ m/^FEATURES/) {
            $feat_flag = 1;
        } elsif ($feat_flag) {
            if ($line =~ m/^\s{1,10}(rRNA)\s+(?:complement\()?<?(\d+)\.\.>?(\d+)/) {
                $feat_flag = $1;
                $rrna_count++;
                $rrna[$rrna_count]->{'start'} = $2;
                $rrna[$rrna_count]->{'end'} = $3;
                # the following filters several Acropora sequences which 
                # include only 9 nucleotides of 18S
                # the rest of the sequence is MT control region
                return 0 if ($rrna[$rrna_count]->{'start'} == 1 &&
                             $rrna[$rrna_count]->{'end'}   == 9);
            } elsif ($line =~ m/^\s{1,10}(CDS|gene|mRNA|ncRNA|CONTIG|misc_RNA|primer_bind|variation)\s+(\S.+)/) {
                $feat_flag = $1;
            } elsif ($line =~ m/^\s{1,10}(\S+)\s+(?:complement\()?<?(\d+)\.\.>?(\d+)/) {
                $feat_flag = $1;
            } elsif ($line =~ m/^\s{20,}(.*)/) {
                if ($feat_flag eq 'rRNA') {
                    $rrna[$rrna_count]->{'misc'} .= $1;
                }
            } else { die "unexpected: $line\n$acc"; }
        } elsif ($line =~ m/^REFERENCE/ && !$hit_ref) {
            $tax_flag = 0;
            $tax_line =~ s/ //g;
            @tax = split /;/, $tax_line;
            $tax_line = '';
            $hit_ref = 1;
        } elsif ($tax_flag) {
            chomp $line;
            $tax_line .= $1 if ($line =~ m/^\s+(\S.+)/);
        } elsif ($seq_flag) {
            if ($line =~ m/^\/\//) {
                $seq_flag = 0;
            } elsif ($line =~ m/^\s+\d+ (.*)/) {
                $seq .= $1;
            }

        } elsif ($line =~ m/^\s+Eukaryota; Metazoa; Cnidaria; (\S+); (\S+);(?: ([^;]+))?/) {
            
            $class = "$1/$2";
            $class .= "/$3" if ($3);
        }
    }
    $seq =~ s/ //g;
#    $seq =~ s/[Nn]//g;
    my $count_ambig = () = $seq =~ m/[^AaCcTtGgUu]/g;
    return 0 if ($count_ambig > $MAX_AMBIG);
    $rh_data->{$org}->{$acc} = {'len' => $len,  'seq'  => $seq, 
                                'tax' => \@tax, 'rrna' => \@rrna};
}

