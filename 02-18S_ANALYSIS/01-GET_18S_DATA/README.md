# 18s Analyses

### 18s.cnidaria.genbank.gz

result of the following GenBank search on 2020, July 8:
```((Cnidaria[ORGN] AND (18S OR "small subunit ribosomal")) BUTNOT Nematostella[ORGN]) OR AF254382```
AF254382 was added.

### get_18S_fasta_from_genbank.pl

script used to extract 18s FASTA from the 18s.cnidaria.genbank.gz file

### filter_based_on_prior_knowledge.pl

script to remove/rename a few sequences

### 18s.filtered.fa

output of filter_based_on_prior_knowledge.pl.  Used to generate 18s tree

### 01-PREFILTERED/18s.fa

output of get_18S_fasta_from_genbank.pl and input to filter_based_on_prior_knowledge.pl

### 01-PREFILTERED/18s.long.fa

output of get_18S_fasta_from_genbank.pl (long deflines)


