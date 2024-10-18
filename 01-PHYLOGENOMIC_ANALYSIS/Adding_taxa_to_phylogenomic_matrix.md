## Adding taxa to our phylogenomic matrix

Below is an untested set of instructions on how to add taxa to our phylogenomic matrix

NOTE: We did not include a partition file in the original repository, but I have added one to this repo (cnido_748_untrimmed_b_partition_info.txt).
Individual alignments can be created using this file and the full concatenated matrix, but to save some time, I copied the pre-concatenated FASTA files (untrimmed) to here:
    http://ryanlab.whitney.ufl.edu/downloads/debiasse_et_al_cnid_partitions_fasta.tar.gz
 
Before doing the following steps. Install hmm2aln.pl from here: https://github.com/josephryan/hmm2aln.pl
It will also tell you to install hmmer and JFR-PerlModules.  Instructions for doing this via conda are provided.
 
To concatenate taxa to this data set, do the following:
 
1.	uncompress/untar the fasta files: 
 
```bash
gzip -dc debiasse_et_al_cnid_partitions_fasta.tar.gz | tar -xvf -
```
 
2.	cd into the newly created directory:
 
```bash
cd debiasse_et_al_cnid_partitions_fasta
```
 
3.	Create hidden Markov models for each alignment:
 
```bash
ls -1 | perl -ne 'chomp; m/(.*).untrimmed.fa.mafft/; print "hmmbuild $1.hmm $_\n";' | sh 
```
 
4.	Use hmm2aln.pl to align your sequences to these newly created hmms:
 
```bash
ls -1 | perl -ne 'chomp; m/(.*).untrimmed.fa.mafft/; print "hmm2aln.pl --hmm=$1.hmm --fasta=YOURFASTAFILE --threads=NUMTHREADS > $1.h2a.fa\n";' | sh
```

Once you have done this, I recommend looking at several of these alignments and making sure your sequences have been aligned in reasonable fashion. I use UGENE to view alignments: https://ugene.net/

NOTE: hmm2aln.pl removes insertions created by your sequence so some trimming may occur in your sequences (unlike the rest of the dataset which has no trimming).  I think this is acceptable, but you should mention it if you end up using this approach.

5. generate a tree for each alignment (we use iqtree - iqtree.org):

```bash
iqtree -s GENEID.h2a.fa
```

6. run phylotreepruner to remove isoforms/paralogs (https://sourceforge.net/projects/phylotreepruner/)

```bash
phylopypruner --output=GENEID.out --dir DIRECTOR_W_TREES_AND_ALNS --mask longest --min-support 0.5 --prune MI
```

7. If you a clear single entry from your new genome great, if not, make sure you have a version of the original alignment (from the tar file I sent) with a FASTA entry for your species with all gaps as sequence

8. Concatenate all the sequences (fasta2phylomatrix is a script that installs when you install my perl modules, which would have been done while installing hmm2aln.pl)

```bash
fasta2phylomatrix --dir=DIR_OF_FASTA_ALIGNMENTS --partition=OUT_PARTITION_FILE > concat_align.fa
```
