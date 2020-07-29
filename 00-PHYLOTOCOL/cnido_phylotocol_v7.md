# PLANNED ANALYSES FOR ASSEMBLING THE CNIDARIAN TREE OF LIFE 
 Principle Investigators: Joseph Ryan and Melissa DeBiasse 
 Draft or Version Number: v.1.1  
 Date: 29 July 2020  
 Note: updates to this document will be tracked through github
 
## 1 INTRODUCTION: BACKGROUND INFORMATION AND SCIENTIFIC RATIONALE  

### 1.1 _Background Information_  
Cnidaria is an exclusively aquatic clade of >13,000 species. Relationships within Cnidaria have been controversial for many years. Here, we outline our planned phylogenetic analyses of an expanded set of cnidarian transcriptomes. 

### 1.2 _Rationale_ 
Three recent phylogenomic studies (Chang et al. 2015; Zapata et al. 2015; Kayal et al. 2018) produced congruent results for most of the main cnidarian lineages. Several lineages, especially within Anthozoa, remain under sampled or not sampled at all. 

### 1.3 _Objectives_   
This project will expand upon previous phylogenomic studies of cnidarians in terms of taxon sampling. We will assemble transcriptomes for several cnidarian taxa, add these data to previously published assemblies, construct a matrix of orthologous loci, and estimate a phylogenetic tree.We will also build a tree of many hundred cnidarian species by performing a phylogenetic analysis of 18S RNA sequences that are scaffolded by our transcriptomic phylogeny. 


## 2 STUDY DESIGN  

#### 2.1 Assemble transcriptomes from species sequenced for this study and from species with data available in the NCBI SRA  

2.1.1 download RNA-Seq data from SRA

```
perl get_sra.pl [ACCESSION_NUMBER]
```

2.1.2 edit deflines in fastq reads

```
perl fix_names.pl 1.fastq 2.fastq unpaired.fq > fix.out 2> fix.err
```

2.1.3 assemble transcriptomes with Trinity v2.8.5

```
Trinity --seqType fq --max_memory [#G] --CPU [#] --trimmomatic --left 1.fastq.renamed --right 2.fastq.renamed --full_cleanup --include_supertranscripts --output cnido_sp.trinity > trin.out 2> trin.err
```

#### 2.2 Translate the nucleotide transcriptomes into amino acid sequences with TransDecoder v3.0.1. We set the –m flag to 50 and used the results from BLAST searches to inform the final TransDecoder prediction step

```
TransDecoder.LongOrfs -t [cnido_sp.trinity.fa] -m 50 > td.out 2> td.err
```

```
blastp -query longest_orfs.pep -db swissprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads [#] > blastp.out 2> blastp.err 
```

```
TransDecoder.Predict -t [cnido_sp.trinity.fa] --retain_blastp_hits out.blastp.out > tdp.out 2> tdp.err
```

#### 2.3 Filter non-target algal symbiont sequences in the program Alien Index (https://github.com/josephryan/alien_index) using a database of metazoan and non-metazoan sequences  

```
blastp -query [cnido_sp.pep.fa] -db ai.fa -outfmt 6 -max_target_seqs 1000 -seg yes -evalue 0.001 -out cnido_sp.blast.out > ai_1.out2> ai_1.err
```

```
./alien_index --blast=cnido_sp.blast.out --alien_pattern=ALIEN [out.alien_index] > ai_2.out 2> ai_2.err 
```

```
remove_aliens.pl [out.alien_index] [cnido_sp.pep.fa] > cnido_sp_no_aliens.pep.fa > ai_3.out 2> ai_3.err
```

2.3.1 rename the definition lines

```
replace_deflines.pl --pad=6 --fasta=cnido_sp_no_aliens.pep.fa --prefix=[Genus_species] > cnido_sp_no_aliens_renamed.pep.fa
```

#### 2.4 Calculate transcriptome completeness with BUSCO v2/v3 implemented by the gVolante webserver (https://gvolante.riken.jp/). We searched each transcriptome against the eukaryote ortholog set  

#### 2.5 Identify orthologs among taxa using OrthoFinder v2.2.3 to identify orthologs among taxa. A custom script to run the diamond searches (ortho_diamond.pl) is available in this repository 

```
orthofinder -f [dir_w_protein_fasta_files] -op > of.out 2> of.err

```

```
diamond blastp -d BlastDB -q species.fa -o species_out.txt -e 0.001 -p [# of cores] -f 6 > species.out 2> species.err &
```

```
orthofinder -b [dir_w_blast_results] -a 16 -M msa -os > ofb.out 2> ofb.err
```

#### 2.6 Generate single copy orthogroups using the script ```filter_ogs_write_scripts.pl``` (available in https://github.com/josephryan/RyanLabPhylogenomicTools). The script allows the user to define the minimum number of taxa and the maximum number of duplicates per taxon allowed per orthogroup (in this study XXX% (XXX taxa) and X duplicates) and after this filtering, automates the following steps:   

2.6.1 sequences within each orthogroup are aligned using Mafft v7.309 

```
mafft-linsi --localpair --maxiterate 1000 --thread 20 [infile] >mafft.out 2> mafft.err
```

2.6.2 alignments are refined using Gblockswrapper v0.03 (https://goo.gl/fDjan6)

```
Gblockswrapper [infile.mafft] > outfile.mafft-gb > gbw.out 2> gbw.err
```

2.6.3 alignments with sequences with more than a set threshold of gaps are removed with the ```remove_n_percent_gaps.pl``` script, available in https://github.com/josephryan/RyanLabPhylogenomicTools. We set the threshold to 50% 

```
remove_n_percent_gaps.pl [outfile.mafft-gb] > rnpg.out 2> rnpg.err
```

2.6.4 maximum-likelihood orthogroup gene trees are estimated in IQTree v1.5.5 

```
iqtree-omp -s [infile.mafft-gb] -nt AUTO -bb 1000 -m LG -pre [output prefix] > iq.out 2> iq.err
```

2.6.5 paralogs are pruned in PhyloTreePruner v1.0

```
java PhyloTreePruner [infile.tree] 28 [infile.align] 0.5 u > ptp.out 2> ptp.err
```

#### 2.7 Concatenate the single-copy loci filtered from step 2.7 to create a matrix and partition file for use in downstream phylogenomic analyses using ```fasta2phylomatrix``` (available in https://github.com/josephryan/RyanLabPhylogenomicTools). Definition lines in each fasta file were edited (```perl -pi.orig -e 's/\|.*$//;' *.fa```) prior to running ```fasta2phylomatrix```  

#### 2.8 Estimate species phylogeny from the concatenated matrix using maximum likelihood 

```
iqtree-omp -s cnido_matrix.fa -pre cnido -spp cnido.nex -nt AUTO -m TEST -bb 1000 > iq.stdout 2> iq.err
``` 

#### 2.9 Estimate the cnidarian species phylogeny from concatenated and untrimmed orthogroup alignments to test the effect of Gblockswrapper in step 2.6.2

2.9.1 modify the deflines in the original orthogroup fasta files 
```
replace_cnid_defs.pl original_orthogroup.fa > original_orthogroup_new_deflines.fa
```

2.9.2 remove sequences that were pruned in step 2.6.5 from the original orthogroup fasta files
```
perl revert_to_untrimmed.pl --original_orthogroup_new_deflines.fa --subset_fasta=original_orthogroup_new_deflines_subset.fa
``` 

2.9.3 align sequences within each pruned orthogroup fasta file
```
mafft original_orthogroup_new_deflines_subset.fa > original_orthogroup_new_deflines_subset.mafft
```

2.9.4 concatenate fasta files generated in step 2.9.3 to create a matrix and partition file using ```fasta2phylomatrix``` (available in https://github.com/josephryan/RyanLabPhylogenomicTools)

2.9.5 estimate ML species tree from untrimmed concatenated alignment 
```
iqtree-omp -s cnido_matrix_untrimmed.fa -pre cnido_untrimmed -spp cnido_untrimmed.nex -nt AUTO -m TEST -bb 1000 > iq.stdout 2> iq.err
```

#### 2.10 Estimate cnidarian phylogeny using 18s sequences to infer species level relationships with higher order relationships constrained to match the topology of the concatenated 748-locus ML phylogeny 

2.10.1 download GenBank records for 18S Cnidaria sequence
```
((Cnidaria[ORGN] AND (18S OR "small subunit ribosomal")) BUTNOT Nematostella[ORGN]) OR AF254382
```

2.10.2 filter records on the following criteria (a) ...
and create a file with FASTA sequences

```
perl get_18S_fasta_from_genbank.pl 18s.cnidaria.genbank > 18s.fa
```

2.10.3 align 18S sequences with ssu-align (Nawrocki and Eddy, 2013).

```ssu-align -f 18s.fa ssu.dir > ssu-align.out 2> ssu-align.err
```

2.10.4 remove positions with low posterior probability of positional homology as calculated by SSU-mask and convert stockholm to fasta (custom script that uses esl-reformat from HMMer)

```
ssu-mask -m /usr/local/ssu-align-0.1.1/env/ssu-align-0p1.cm ssu.dir > ssu-mask.out 2> ssu-mask.err
stockholm2fasta.pl ssu.dir/ssu.dir.eukarya.mask.stk > ssu.dir/ssu.dir.eukarya.mask.stk.fa
```

2.10.5 construct constraint tree from transcriptome-based phylogeny, by renaming taxa to match 18s names and removing taxa that are not in both datasets.

```
construct_constraint_tree.pl transcriptome.tre > transcriptome.constraint.tre
```

2.10.6 unroot the constraint tree (I think iqtree balks at constraint tree unless tree is specifically unrooted; if not, than will skip)

```
unroot.R
```

2.10.7 construct a phylogenetic tree 

```
iqtree -nt AUTO -s ssu.dir/ssu.dir.eukarya.mask.stk.fa -g transcriptomic_constraint.v4.unrooted.tre -m TEST > iq.out 2> iq.err
```

## 3 WORK COMPLETED SO FAR

13 July 2020: 2.1-2.8, 2.9.1-2.9.4 

## 3 PROGRAMS REFERENCED  

Emms, D. M., & Kelly, S. (2015). OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biology, 16(1), 157. 

Gblockswrapper: http://bit.ly/2svaKcR

Kocot, K. M., Citarella, M. R., Moroz, L. L., & Halanych, K. M. (2013). PhyloTreePruner: a phylogenetic tree-based approach for selection of orthologous sequences for phylogenomics. Evolutionary Bioinformatics Online, 9, 429.

Lartillot, N., Rodrigue, N., Stubbs, D., & Richer, J. (2013). PhyloBayes MPI: phylogenetic reconstruction with infinite mixtures of profiles in a parallel environment. Systematic Biology, 62(4), 611-615.

Nawrocki EP, Eddy SR. Infernal 1.1: 100-fold faster RNA homology searches. Bioinformatics. 2013 Nov 15;29(22):2933-5.

Nishimura O, Hara Y, Kuraku S. (2017) gVolante for standardizing completeness assessment of genome and transcriptome assemblies. Bioinformatics. 33(22), 3635-3637.

Nguyen, L. T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2014). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. Molecular Biology and Evolution, 32(1), 268-274.

TransDecoder: https://transdecoder.github.io/

Yamada, K. D., Tomii, K., & Katoh, K. (2016). Application of the MAFFT sequence alignment program to large data—reexamination of the usefulness of chained guide trees. Bioinformatics, 32(21), 3246-3251.

## APPENDIX

Version : Date : Significant Revisions  
1.1 29 July 2020 Before identifying single-copy orthogroups (step 2.6), we removed five cnidarian species with a high number of duplicates per core gene (Anemonia viridis, Heliopora coerulea, Montastrea faveolata, Periphylla periphylla, Stomphia coccinea). In the first ML tree, Heteractis crispa grouped with the avian outgroup and Muricea muricata grouped with Eunicea due to contamination. Further investigation revealed contamination in these taxa so we removed them from the analysis and re-estimated the ML tree.   
1.2  
1.3  
1.4 
