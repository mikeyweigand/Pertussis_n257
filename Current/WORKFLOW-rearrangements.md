## Workflow for rearrangement analysis:
---
1. Group genomes into collinear clusters.
 + Exhaustive pairwise alignment with __progressiveMauve__ ([Darling *et al.* 2010](http://www.ncbi.nlm.nih.gov/pubmed/20593022)).
    + __pairwise-mauve-all.pl__
    + __pairwise-mauve-parallel.sh__
 + Identify collinear genome pairs.
    + __mauve.backbone-clean.pl__ (remove IS-containing gaps)
    + __mauve.collinear-check.pl__
    + Concatenate output summaries.
 + Cluster genomes using __mcl__ ([Markov Cluster Algorithm](http://http://micans.org/mcl/)).
	+ __mauve.collinear-mcl.pl__
1. Align non-redundant subset of structures.
 + Exctract subset:  `cut -f3 mcl-clusters.txt`
 + __progressiveMauve__ with optimized parameters.
 + __Scriptbox/src/mauve.backbone2matrix.pl__
1. Map permutation matrices to full dataset of collinear genomes.
 + __Scriptbox/src/mauve.matrix-transfer.pl__
1. Cluster genomes with __MLGO__ ([Hu *et al.* 2014](https://www.ncbi.nlm.nih.gov/pubmed/25376663)).
 + View, annotate tree with [iTOL](http://itol.embl.de/).
