## Workflow for matching variable IS481 insertions:  
For matching variable IS481 insertions in a set of genomes as described in __Weigand *et al.* 201X__  

1. Exhaustive pairwise alignment of genomes with progressiveMauve ([Darling *et al.* 2010](http://www.ncbi.nlm.nih.gov/pubmed/20593022)).
1. [BLASTn](http://www.ncbi.nlm.nih.gov/books/NBK279671/) search of IS481 ([M22031](http://www.ncbi.nlm.nih.gov/nuccore/144060/)) insertions in each genome.  
  * Use: -outfmt '6 qseqid sseqid qlen length qstart qend sstart send'
1. Match IS481 insertions between genomes in each pairwise alignment.  
  * __mauve.backbone-IStracker.pl__
1. Convert IS matrices to rbm format.
  * __ISmatrix2rbm.pl__  
1. Combine matches across all genomes.  
  * __ogs.rb__ ([enveomics](http://enve-omics.ce.gatech.edu/enveomics/)).
1. Identify and add unmatched, strain-specific insertions.
  * __ISmatrix.add-strain-spec.pl__
1. Final consistency check for unidentified duplications.
  * __ISmatrix.duplication-check.pl__


__WARNING 1:__ Many scripts in this workflow are currently written to specifically recognize name formats containing PDL strain-ids and GenBank accession numbers, eg. 'A001_CP000001'.

__WARNING 2:__ In principle, this approach should work for any set of complete genome assemblies. However, it is limited by the accuracy of progressiveMauve and, so in practice, suffers greatly as alignment complexity increases due to multiple rearrangements. Without an upstream method for alignment correction, this approach is best used only for analyzing sets of collinear genomes. It may also be suitable for datasets with only a few inversions, but has not been extensively tested.  
