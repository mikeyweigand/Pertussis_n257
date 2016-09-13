#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

#description: create signed permutation matrix from progressiveMauve *.backbone file.

&GetOptions(	'backbone=s' => \my$bbone,	# .backbone file
		'xmfa=s' => \my$xmfa,		# .xmfa
		'num=s' => \my$mingen,		# minimum number of genomes required for LCB
		'minlen=s' => \my$minlen,	# minimum LCB length
		'out=s' => \my$out);		# name prefix for output table/matrix

($bbone and $xmfa and $mingen and $minlen and $out) or &HELP_MESSAGE;


#store ordered list of genomes from xmfa file
my@genomes=();
open XMFA, "$xmfa";
while(my$x = <XMFA>){
	chomp$x;
	if($x =~ /\#Sequence[0-9]+File/){
		my@sx = split("\t",$x);
		my$g = basename($sx[1]);
		$g =~ s/\.[a-z]+$//;
		push(@genomes,$g);

	}
}
close XMFA;

#open .backbone file and store homologous block coordinates
my%hash=();
my%ghash=();
my%bhash=();
my$bnum = 1;

open BBONE, "$bbone";
while(my$bb = <BBONE>){
	chomp$bb;
	my@sbb = split("\t",$bb);
	my$gnum = 0;
	unless($bb =~ /^seq/){

		$hash{ $bnum } = 0;

		for (my$ii=0; $ii < @sbb; $ii+=2){			## only pass blocks >minlen
			if(abs($sbb[$ii]-$sbb[$ii+1]) >= $minlen){	
				
				my$blen = abs($sbb[$ii]-$sbb[$ii+1]);
				
				#count number of genomes containing block >= minlen
				$hash{$bnum}++;
				#store block coordinate information for each genome
				@{ $ghash{ $genomes[$ii/2] }{abs($sbb[$ii])}} = ($bnum,$sbb[$ii],$sbb[$ii+1],$blen);
				@{ $bhash{ $bnum }{ $genomes[$ii/2] } } = ($sbb[$ii],$sbb[$ii+1],$blen);
				#print $bnum."--".$genomes[$ii/2]."\n";
			}
		}
		$bnum++;
	}
	
}
close BBONE;

#exclude blocks present in less than minimum number of genomes
foreach my$j (sort{$a <=> $b}(keys%hash)){
	unless($hash{$j} >= $mingen){	
		delete($hash{$j});
	}else{
	#print $j."\n";
	}
}


#create signed permutation matrix for MLGO
my$matrix=$out.".blocks-matrix.txt";
open MATRIX, ">$matrix";
my@lengths=();

foreach my$g2 (@genomes){
	print MATRIX ">".$g2."\n";
	my$glen=0;
	foreach my$cc (sort{$a <=> $b}keys$ghash{$g2}){
		if(exists($hash{ $ghash{$g2}{$cc}[0] } )){
			
			#if coordinates are negative, block is inverted (negative sign)
			if( $ghash{$g2}{$cc}[1] >= 1){
				print MATRIX $ghash{$g2}{$cc}[0] . " ";
			}else{
				print MATRIX $ghash{$g2}{$cc}[0]*-1 . " ";
			}
			$glen+=$ghash{$g2}{$cc}[-1];
		}
	}
	print MATRIX "\$\n";
	push(@lengths, $glen);
}
close MATRIX;

my@slengths = sort{$a <=> $b}@lengths;

#create table of block coordinates and lengths
my$table = $out.".blocks-table.txt";
open TABLE, ">$table";
print TABLE "block";
foreach my$h (@genomes){
	print TABLE "\t".$h."_left\t".$h."_right\t".$h."_len";
}
print TABLE "\n";

foreach my$t (sort{$a <=> $b}(keys%hash)){
	print TABLE "block_".sprintf("%05d",$t);
	#print $t."\n";
	foreach my$g3 (@genomes){
		#print "\t".$g3."\n";
		if( exists( $bhash{ $t }{ $g3 } )){
			print TABLE "\t". join("\t", @{ $bhash{ $t }{ $g3 } });
		}else{
			print TABLE "\t0\t0\t0";
		}
	}
	print TABLE "\n";
}


#output summary of results
print "\n\tThere are ".scalar(@genomes)." genomes with ".scalar(keys%hash)." homologous blocks that meet the following cut-offs:\n";
print "\t\tPresent in at least $mingen genomes\n";
print "\t\tWith length >= $minlen bp\n";
print "\t\tWhich comprise at least $slengths[0] bp per genome\n\n";

	
sub HELP_MESSAGE { die "
.Description:
   Extracts block coordinates from a progressiveMauve .backbone and .xmfa output file and generates:
   	(1) coordinate table, with block lengths
   	(2) signed permutation matrix for input to MLGO.
   
.Usage: $0 -backbone <*.backbone> -xmfa <*.xmfa> -num 5 -minlen 500 -out <prefix>
   
   [mandatory options]
   -backbone	Input backbone alignment file.
   -xmfa	Input xmfa alignment file.
   -num		Minimum number of genomes required per LCB (example: 5).
   -minlen	Minimum LCB length in bp (example: 500).
   -out		Filename prefix for output.

" }

