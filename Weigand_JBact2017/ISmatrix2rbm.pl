#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;

&GetOptions(	'indir=s' => \my$indir,		# input directory of matrices
		'rbmdir=s' => \my$rbmdir,		# output directory of rbm files
		'nondir=s' => \my$nondir);		# output directory of unmatched insertions

($indir and $rbmdir and $nondir) or &HELP_MESSAGE;

unless(-d $rbmdir){
	system("mkdir $rbmdir");
}else{
	print "\nERROR: '$rbmdir' already exists. Please define a new directory for -rbmdir.\n";
	&HELP_MESSAGE;
}
unless(-d $nondir){
	system("mkdir $nondir");
}else{
	print "\nERROR: '$nondir' already exists. Please define a new directory for -nondir.\n";
	&HELP_MESSAGE;
}


while(my$i = <"$indir/*matrix.txt">){
	print $i."\n";

	open IN, "$i";
	my$rbm = basename($i, "-matrix.txt");
	my@ids = split("-vs-",$rbm);
	$rbm =~ s/-vs//;
	#print $rbm."\t".join("--",@ids)."\n";
	open RBM, ">$rbmdir/$rbm.rbm";
	open NON, ">$nondir/$rbm-unmatched.txt";

	while(my$j = <IN>){
		unless($j =~ /^IS/){
			chomp$j;
			my@sj=split("\t",$j);

			unless(($sj[4] =~ /NA/) || ($sj[5] =~ /NA/)){
				my($ise) = $sj[0] =~ /(?<=\d_)(IS\d+)/;
				#print $ise."\n";
				print RBM $ids[0].".".$ise.".".$sj[4]."\t".$ids[1].".".$ise.".".$sj[5]."\t100.0\n";
			}else{
				my($ise) = $sj[0] =~ /(?<=\d_)(IS\d+)/;
				if(($sj[4] =~ /NA/)){
					print NON $ids[1].".".$ise.".".$sj[5]."\n"
				}
				if(($sj[5] =~ /NA/)){
					print NON $ids[0].".".$ise.".".$sj[4]."\n"
				}
			}
		}

	}
	close RBM;
	close NON;
}

##########################
sub HELP_MESSAGE { die "
.Description:
   Takes a directory of pairwise-matched IS matrices, such as those output from 'mauve.backbone-IStracker.pl', and converts each to *.rbm files suitable for combining with 'ogs.rb'.

.Usage: $0 -indir -rbmdir -nondir

   [mandatory]
   -indir 	<path>	Input directory of matrices from 'mauve.backbone-IStracker.pl'.
	 		Expects files within to be named: 'genome1-vs-genome2-matrix.txt'
   -rbmdir	<path>	New directory of output rbm files.
	 		Will be named like: 'genome1-genome2.rbm'
   -nondir	<path>	New directory of unmatched IS insertions.
	 		Will be named like: 'genome1-genome2-unmatched.txt'

" }
