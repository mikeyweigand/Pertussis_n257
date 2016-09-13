#!/usr/bin/perl -w
use strict;
use Getopt::Long;

&GetOptions(	'in=s' => \my$matrix,		# input permutation matrix
		'dir=s' => \my$dir,		# output directory
		'out=s' => \my$out);		# output file adjacency table

($matrix and $dir and $out) or &HELP_MESSAGE;

system("mauve.MLBE-adjacencies.pl $matrix $dir");
system("Table.merge.pl $dir/*adj.txt > $dir/tmp.txt");


my%ahash=();
$ahash{ -1 } = "\$";

open BLOCKS, "$dir/block-coords.txt";
while(my$b = <BLOCKS>){
	chomp$b;
	my@sb=split("\t",$b);
	my@coords=split("-",$sb[1]);
	$ahash{ $coords[0] } = $sb[0];
	$ahash{ $coords[1] } = $sb[0];
}
close BLOCKS;
open OUT, ">$out";
open ADJ, "$dir/tmp.txt";
while(my$adj = <ADJ>){
	if($adj =~ m/^Tag/){
		$adj =~ s/_adj//g;
		print OUT $adj;
	}else{
		chomp $adj;
		my@sadj=split("\t",$adj);
		(my$cc = $sadj[0]) =~ s/Adjacency_//;
		my@scc=split(",", $cc);
		
		print OUT $sadj[0]." [LCBs_".$ahash{ $scc[0] }.",".$ahash{ $scc[1] }."]\t";
		shift@sadj;
		print OUT join("\t",@sadj)."\n";
	}

}
close OUT;
close ADJ;

sub HELP_MESSAGE { die "
.Description:	
   Takes an MLGO-formated permutation matrix and outputs a table of specific 'adjaceny' presence/absence in each genome. This script is simply a 'wrapper' and depends on 'mauve.MLBE-adjacencies.pl' and 'Table.merge.pl' to do the real work, both should be in your PATH.

.Usage: $0 -in in.txt -dir ./path -out out.txt
   
   [mandatory options]
   -in		Input MLGO-format permutation matrix, like the output from 'mauve.backbone2matrix.pl'
   -dir		Output directory to store individual adjacency lists.
   -out		Filename for output table.
   
   [dependencies]
   mauve.MLBE-adjacencies.pl	(adapted from MLGO; http://www.geneorder.org/)
   Table.merge.pl		(from the Enveomics collection; http://enve-omics.ce.gatech.edu/enveomics/)

" }
