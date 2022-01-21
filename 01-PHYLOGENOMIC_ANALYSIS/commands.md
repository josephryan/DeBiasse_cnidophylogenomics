### transcriptome assembly with Trinity

time trinityrnaseq-Trinity-v2.8.5/Trinity --seqType fq --max_memory 750G --CPU 45 --trimmomatic --left R1.fq.gz --right R2.fq.gz --full_cleanup --include_supertranscripts --output transcriptome > trin.out 2> trin.err

### translate transcriptome into protein sequence with TransDecoder

TransDecoder-3.0.1/TransDecoder.LongOrfs -t transcriptome.Trinity.SuperTrans.fasta
blastp -query transcriptome.Trinity.SuperTrans.fasta.transdecoder_dir/longest_orfs.pep -db swissprot -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 39 > blastp.out 2> blast.err
TransDecoder-3.0.1/TransDecoder.Predict -t transcriptome.Trinity.SuperTrans.fasta --retain_blastp_hits Trinity.SuperTrans.fasta.transdecoder_dir.blastp.out > predict.out 2> predict.err

### filter non-target, non-metazoan (i.e. contaminant) sequences with Alien Index

diamond blastp --query transcriptome.Trinity.SuperTrans.fasta.transdecoder.pep.fa --db ai_db_diamond --outfmt 6 --max-target-seqs 1000 --seq yes --evalue 0.001 --threads 6 --out ai_diamond.out > diamond.stdout 2> diamond.err
alien_index --ai_diamond.out --alien_pattern=ALIEN_ > ai.out 2> ai.err
perl remove_aliens.pl transcriptome.Trinity.SuperTrans.fasta.transdecoder.pep.fa > Trinity.SuperTrans_no_aliens.pep.fa 2> ai2.err
replace_deflines.pl --pad=6 --fasta=transcriptome.Trinity.SuperTrans_no_aliens.pep.fa --prefix=[Genus_species] > transcriptome.Trinity.SuperTrans_no_aliens.renamed.pep.fa

### identify orthology groups across transcriptomes with OrthoFinder

orthofinder -f /dir_with_transcriptome.Trinity.SuperTrans_no_aliens.renamed.pep.fa_files -op > blast_commands.txt
diamond blastp -d BlastDB -q transcriptome.Trinity.SuperTrans_no_aliens.renamed.pep.fa -o transcriptome.Trinity.SuperTrans_no_aliens.renamed.pep.fa_out.txt -e 0.001 -p [# of cores] -f 6 > ortho_blastp.out 2> ortho_blastp.err
orthofinder -b /dir_w_blast_results -a 16 -M msa -os > ortho.out 2> ortho.err

### filter orthogroups according to user-defined parameters to retain single-copy orthologs

filter_ogs_write_scripts --ogdir=DIR_W_OG_FA_FILES --name=NAME4OUTFILES --min_sp=MIN_NUM_SPECIES --num_scripts=NUM_SCRIPTS --max_sp_occur=MAXNUMFOREACHSPECIES --threads=NUMTHREADS

### concatenate single-copy orthologs into a matrix and generate partition file

perl remove_numbers_before_fast2phylomatrix.pl INDIR OUTDIR
fasta2phylomatrix  {--fasta=FILE1 --fasta=FILE2 ... or --dir=DIR_OF_FASTAFILES} --partition=OUT_PARTITION_FILE > concat_align.fa

### estimate maximum likelihood tree from concatenated matrix with IQTREE

iqtree-omp -s concat_align.fa -pre concat_align -spp concat_align.nex -nt AUTO -m TEST -bb 1000 > iq.stdout 2> iq.err
