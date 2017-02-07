## Workflow for rearrangement analysis:  
For characterizing genome rearrangements in complete assemblies as described in __Weigand *et al.* 201X__  

1. Multiple genome alignment with __progressiveMauve__ ([Darling *et al.* 2010](http://www.ncbi.nlm.nih.gov/pubmed/20593022))  
  * Align a non-redundant subset of genomes, representing each unique architecture (predetermined).  
  * Parameters used in Weigand *et al.* : --seed-weight=16 --hmm-identity=0.85
2. Calculate permutation matrix from progressiveMauve output.  
  * __mauve.backbone2matrix.pl__
  * Parameters used in Weigand *et al.* : -num 1 -minlen 1500
3. Map permutation matrices to collinear genomes in full dataset.  
  * __mauve.matrix-transfer.pl__
4. Run MLGO ([Hu *et al.* 2014](http://www.ncbi.nlm.nih.gov/pubmed/25376663))  
  * View results in your favorite tree application or with [iTOL](http://itol.embl.de).  
5. Convert permutation matrix (from step 3) into presence/absence table of specific structural 'adjacencies'.  
  * __mauve.matrix2adj.pl__
6. Compare the presence of 'adjacencies' between two lists of genomes and identify those significantly associated with one group vs the other.  
  * __mauve.adj-filter.fisher.pl__
