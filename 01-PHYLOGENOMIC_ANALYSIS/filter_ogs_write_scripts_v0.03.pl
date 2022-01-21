#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# Copyright (C) 2017,2018; Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# for more info see <http://www.gnu.org/licenses/>

use strict;
use warnings;
use JFR::Fasta;
use Pod::Usage;
use Getopt::Long;
use IO::File;
use Data::Dumper;
use POSIX qw/strftime/;

our $VERSION = 0.03;
our $DELIM   = 'ZZZ'; # temp delim should not be in deflines
our $REMOVE_N_PERCENT_GAPS = 'remove_n_percent_gaps.pl';

MAIN: {
    my $rh_opts = get_opts();
    my $min_sp = $rh_opts->{'min_sp'};
    my $max_sp_occur = $rh_opts->{'max_sp_occur'} || 0;
    my $dir = $rh_opts->{'ogdir'};
    my $name = $rh_opts->{'name'};
    my $threads = $rh_opts->{'threads'} || 1;
    (-d $dir) or die "$rh_opts->{'ogdir'} is not a readable directory\n";
    opendir DIR, $dir, or die "cannot open $dir:$!";
    my @seq_files = grep { /\.fa$/ } readdir DIR;
    my $timestr = strftime('%Y-%m-%d.%H-%M-%S',localtime);
    my $script_name = "$name.script.$timestr";
    my $outdir = "$name.$timestr";
    my $stderr_dir = $outdir . '.stderr';
    die "$outdir exists" if (-d $outdir);
    mkdir $outdir or die "cannot mkdir $outdir:$!";
    die "$stderr_dir exists" if (-d $stderr_dir);
    mkdir $stderr_dir or die "cannot mkdir $stderr_dir:$!";
    opendir DIR, $dir or die "cannot opendir $dir:$!";
    my @data = ();
    FA: foreach my $sf (@seq_files) {
        $sf =~ m/(.*).fa/ or die "unexpected: $sf";
        my $og = $1;
        my $count = 0;
        my %species = ();
        my $fp = JFR::Fasta->new("$dir/$og.fa");
        while (my $rec = $fp->get_record()) {
            $count++;
            $rec->{'def'} =~ m/^>([^_]+_[^_]+)_/;
            my $sp = $1;
            $species{$sp}++;
        }
        my $sp_count = scalar(keys %species);
        next unless ($sp_count >= $min_sp);
        next if ($max_sp_occur && exceed_max_sp_occur(\%species,$max_sp_occur));
        my $ratio = $count / $sp_count;
        push @data, [$og,$count,$sp_count,$ratio];
    }
    my $ra_fh = get_handles($rh_opts->{'num_scripts'},$script_name);
    my $fh_count = 0;
    foreach my $ra_d (sort { $a->[3] <=> $b->[3] } @data) {
#next unless ($ra_d);
#print Dumper $ra_d;
        my $og = $ra_d->[0];
        open IN, "$dir/$og.fa" or die "cannot open $dir/$og.fa$!";
        open OUT, ">$outdir/$og.fa" or die "cannot open >$outdir/$og.fa:$!";
        while (my $line = <IN>) {
            print OUT $line;
        }
        close OUT;
        close IN;
        my $curfh = $ra_fh->[$fh_count];
        print $curfh "echo $ra_d->[0] >> $stderr_dir/gbwrapper.out\n";
        print $curfh "echo $ra_d->[0] >> $stderr_dir/gbwrapper.err\n";
        print $curfh "echo $ra_d->[0] >> $stderr_dir/remove_n_percent_gaps.err\n";
        print $curfh "echo $ra_d->[0] >> $stderr_dir/mafft.err\n";
        print $curfh "echo $ra_d->[0] >> $stderr_dir/IQ.out\n";
        print $curfh "echo $ra_d->[0] >> $stderr_dir/IQ.err\n";
        print $curfh "echo $ra_d->[0] >> $stderr_dir/ptp.out\n";
        print $curfh "echo $ra_d->[0] >> $stderr_dir/ptp.err\n";

        print $curfh "mafft --localpair --maxiterate 1000 --thread $threads $outdir/$og.fa > $outdir/$og.mafft 2>> $stderr_dir/mafft.err\n";
        print $curfh "Gblockswrapper $outdir/$og.mafft >> $stderr_dir/gbwrapper.out 2>>$stderr_dir/gbwrapper.err\n";
        print $curfh q~perl -pi -e 's/ //g' ~ . "$outdir/$og.mafft-gb\n";
        print $curfh q~perl -pi -e 's/C_rubr_no_aliens_pep_Gene.(\d+).*/C_rubr.$1/' ~ . "$outdir/$og.mafft-gb\n";
        print $curfh q~perl -pi -e 's/^>.*_pep_(.*)/>$1/' ~ . "$outdir/$og.mafft-gb\n";
        print $curfh "perl -pi -e 's/\\\./$DELIM/' $outdir/$og.mafft-gb\n";

        my $mpg = $rh_opts->{'min_percent_gaps'};
        print $curfh "$REMOVE_N_PERCENT_GAPS $outdir/$og.mafft-gb $mpg > $outdir/$og.mafft-gb.${mpg}_percent_gaps 2>>$stderr_dir/remove_n_percent_gaps.err\n";
        my $og_outdir = "$outdir/$og.IQ_out";
        print $curfh "mkdir $og_outdir\n";
        print $curfh "iqtree-omp -s $outdir/$og.mafft-gb.${mpg}_percent_gaps -nt AUTO -bb 1000 -m LG -pre $og_outdir/$og.iq >> $stderr_dir/IQ.out 2>> $stderr_dir/IQ.err\n";
        print $curfh "perl -pi -e 's/$DELIM/\|/g' $og_outdir/$og.iq.contree $outdir/$og.mafft-gb.${mpg}_percent_gaps\n";
#        print_substitutions($curfh);
        print $curfh "java PhyloTreePruner $og_outdir/$og.iq.contree $min_sp $outdir/$og.mafft-gb.${mpg}_percent_gaps 0.5 u >> $stderr_dir/ptp.out 2>> $stderr_dir/ptp.err\n\n";

        if ($fh_count == 0) {
            $fh_count++;
        } elsif (($fh_count + 1) % $rh_opts->{'num_scripts'}) {
            $fh_count++;
        } else {
            $fh_count = 0;
        }
    }
}
#sub print_substitutions {
#    my $curfh = shift;
#
#    print $curfh q~perl -pi -e 's/([A-Z][a-z]{3}_[a-z0-9]+)_([0-9]+\.[0-9]+)/$1\|$2/g' ~;
#    print $curfh "*.IQ_out/*.iq.contree\n";
#    print $curfh q~perl -pi -e 's/([^_]+_[^_]+)_([^_]+\.[^_]+)/$1\|$2/g' ~;
#    print $curfh "$outdir/*.mafft-gb.res\n";
#}

sub exceed_max_sp_occur {
    my $rh_dat = shift;
    my $max    = shift;
    foreach my $count (values %{$rh_dat}) {
        return 1 if ($count > $max);
    }
    return 0;
}

sub get_handles {
    my $num = shift;
    my $script_name = shift;
    my @handles = ();
    for (my $i = 0; $i < $num; $i++) {
        my $name = "$script_name.$i";
        my $fh = IO::File->new($name,'w');
        push @handles, $fh;
    }
    return \@handles;
}

sub get_opts {
    my $rh_opts = {'num_scripts' => 0 };
    my $res = GetOptions ("num_scripts=i" => \$rh_opts->{'num_scripts'},
                     "min_percent_gaps=i" => \$rh_opts->{'min_percent_gaps'},
                               "min_sp=i" => \$rh_opts->{'min_sp'},
                         "max_sp_occur=i" => \$rh_opts->{'max_sp_occur'},
                                "ogdir=s" => \$rh_opts->{'ogdir'},
                                 "name=s" => \$rh_opts->{'name'},
                                   "help" => \$rh_opts->{'help'},
                               "version"  => \$rh_opts->{'version'},
                              "threads=i" => \$rh_opts->{'threads'} );

    die "filter_ogs_write_scripts version $VERSION\n" if ($rh_opts->{'version'});
    pod2usage({-exitval => 0, -verbose => 2}) if($rh_opts->{'help'});
    print "missing --min_sp\n" unless ($rh_opts->{'min_sp'});
    usage() unless ($rh_opts->{'min_sp'});
    print "missing --min_percent_gaps\n" unless ($rh_opts->{'min_percent_gaps'});
    usage() unless ($rh_opts->{'min_percent_gaps'});
    print "missing --num_scripts\n" unless ($rh_opts->{'num_scripts'});
    usage() unless ($rh_opts->{'num_scripts'});
    print "missing --ogdir\n"  unless ($rh_opts->{'ogdir'});
    usage() unless ($rh_opts->{'ogdir'});
    print "missing --name\n"  unless ($rh_opts->{'name'});
    usage() unless ($rh_opts->{'name'});

    return $rh_opts;
}

sub usage {
    die "usage: $0 --ogdir=DIR_W_OG_FA_FILES --name=NAME4OUTFILES --min_sp=MIN_NUM_SPECIES --num_scripts=NUM_SCRIPTS --min_percent_gaps=MIN_PERCENT_GAPS [--max_sp_occur=MAXNUMFOREACHSPECIES] [--threads=NUMTHREADS] [--version]\n";
}

__END__

=head1 NAME

B<filter_ogs_write_scripts> - write scripts to filter orthofinder outputs

=head1 AUTHOR

Joseph F. Ryan <joseph.ryan@whitney.ufl.edu>

=head1 SYNOPSIS

filter_ogs_write_scripts --ogdir=DIR_W_OG_FA_FILES --name=NAME4OUTFILES --min_sp=MIN_NUM_SPECIES --num_scripts=NUM_SCRIPTS [--max_sp_occur=MAXNUMFOREACHSPECIES] [--threads=NUMTHREADS] [--version]

=head1 OPTIONS

=over 2

=item B<--ogdir>

directory with orthogroup fasta files (after running orthofinder)

=item B<--name>

name used for output directories and scripts

=item B<--min_sp>

minimum number of species in an orthogroup. Orthogroups not meeting this requirement will not show up in the outdir or scripts.

=item B<--num_scripts>

divide jobs up into this many scripts. Allows chunks to be spread across servers

=item B<--max_sp_occur=MAXNUMFOREACHSPECIES>

maximum number of occurrences of any species in a particular orthogroup. Orthogroups that exceed this will not show up in the outdir or scripts.

=item B<--help>

Print this manual

=item B<--version>

Print the version. Overrides all other options.

=back

=head1 DESCRIPTION

Will create scripts that can be used to process orthogroups towards building a phylogenomic matrix

=head1 BUGS

Please report them to the author

=head1 COPYRIGHT

Copyright (C) 2017,2018; Joseph F. Ryan

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
