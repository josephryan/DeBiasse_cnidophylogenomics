# Make Constraint Tree

### Generate initial constraint tree by matching phylogenomics and 18s data

```perl construct_constraint_tree.pl > constraint.v3.tre```

### Produce unrooted tree (constraint.v3.unrooted.tre)

```Rscript unroot.R constraint.v3.tre```

### Make Cubozoa a polytomy (to mirror conflict in our phylogenomic trees)

```perl -pi -e 's/\(Alatina_alata_CU,Chironex_fleckeri_CU\)/Alatina_alata_CU,Chironex_fleckeri_CU/' constraint.v3.unrooted.tre```


