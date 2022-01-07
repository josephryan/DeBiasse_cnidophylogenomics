#!/usr/local/bin/Rscript

library(ape)

# we temporarily root to make sure that the Cubozoa clade
#   is not split so that it can be easly collapsed into a polytomy
tmproot <- "Montastraea_cavernosa_HE"

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("usage: unroot.R NEWICK_FILE", call. = FALSE)
}

tr<-ape::read.tree(args[1]); 
trr<-root(tr,tmproot)

#is.rooted(tr)
unroot<-unroot(trr)
#is.rooted(unroot)
write.tree(unroot, "constraint.v3.unrooted.tre")

