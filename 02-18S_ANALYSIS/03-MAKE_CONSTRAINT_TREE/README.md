# Make Constraint Tree

### commands.md

commands for repeating this analysis

### construct_constraint_tree.pl

script to match 18S and phylogenomic names and produce initial constraint tree
REQUIRES: Phyutility: https://github.com/blackrim/phyutility
REQUIRES: JFR-PerlModules: https://github.com/josephryan/JFR-PerlModules

```perl construct_constraint_tree.pl```

### unroot.R

script to unroot the initial constraint tree (important for downstream steps)

REQUIRES: R -- https://www.r-project.org
REQUIRES: ape (R package): `R -e "install.packages('ape')"`

```Rscript unroot.R```

### cnidophylo_names_abbrev.csv

data file used by construct_constraint_tree.pl

### constraint.v3.tre

output of construct_constraint_tree.pl

### constraint.v3.unrooted.tre

output of unroot.R and perl -pi (see commands.md; constraint for 18s phylogeny)


