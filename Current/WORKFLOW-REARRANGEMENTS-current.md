## Workflow for rearrangement analysis:
---
1. Group genomes into collinear clusters.
 + Exhaustive pairwise alignment with __progressiveMauve__ ([Darling *et al.* 2010](http://www.ncbi.nlm.nih.gov/pubmed/20593022)).
    + __pairwise-mauve-all.pl__
    + __pairwise-mauve-parallel.sh__
 + Identify collinear genome pairs.
    + __mauve.backbone-clean.pl__ (remove IS-containing gaps)
    + __mauve.collinear-check.pl__
    + Concatenate output summaries: e.g. `cat *.sum > all.sum`
 + Cluster genomes using __mcl__ ([Markov Cluster Algorithm](http://http://micans.org/mcl/)).
	+ __mauve.collinear-mcl.pl__
1. Align non-redundant subset of structures.
 + Exctract subset:  `cut -f3 mcl-clusters.txt`
 + __progressiveMauve__ with your optimized parameters.
 + __mauve.backbone2matrix.pl__
1. Map permutation matrices to full dataset of collinear genomes.
 + __mauve.matrix-transfer.pl__
1. Cluster genomes into a tree with __MLGO__ ([Hu *et al.* 2014](https://www.ncbi.nlm.nih.gov/pubmed/25376663)).
 + View, annotate tree with [iTOL](http://itol.embl.de/).

---
NOTES and WARNINGS:
+ The expected input for this workflow is a collection of complete (1-contig) bacterial genome assemblies, each saved in a separate fasta file. Intermediate steps use filename parsing "GenomeA-GenomeB" and so the original fasta filenames must contain only [A-Za-z0-9_].
+ To date, these scripts have only been tested with *Bordetella* species and may require modification for use with other species.
