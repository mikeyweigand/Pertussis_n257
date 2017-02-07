#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Sort::Key::Natural qw(natsort);
use Statistics::R;
use File::Basename;

## default parameters ##
my$p = 0.05;

&GetOptions(	'in=s' => \my$in,		# input matrix
		'l1=s' => \my$list1,		# list 1 of genomes
		'l2=s' => \my$list2,		# list 2 of genomes
		'p=s' => \$p,
		'c' => \my$cor,
		'q' => \my$q,
		'stats=s' => \my$stats,		# output table of summary statistics
		'out=s' => \my$out);		# output matrix

($in and $list1 and $list2 and $out and $stats) or &HELP_MESSAGE;

## store genome lists ##
my%genomes1=();
open LIST, "$list1";
while(my$l = <LIST>){
	chomp$l;
	$genomes1{ $l } = 1;
}
close LIST;

my%genomes2=();
open LIST2, "$list2";
while(my$l2 = <LIST2>){
	chomp$l2;
	$genomes2{ $l2 } = 1;
}
close LIST2;

unless($q){
	print "  List #1 genomes = ". (keys%genomes1)."\n";
	print "  List #2 genomes = ". (keys%genomes2)."\n";
	print "  Matched indices:\n";
}

## filter input matrix ##

open IN, "$in";
open OUT, ">$out";
my@ind1=(); #list1
my@ind2=(); #list2
my@listnames=('List');
my@others=(); 
my$c=0;
my$j=0;
my%adjhash=();
my%adjrows=();
my$R = Statistics::R->new();
while(my$i = <IN>){
	## find indices for listed genomes in matrix header ##
	if($c == 0){
		print OUT $i; #print header to output. 
		chomp$i;
		my@header=split("\t",$i);
		for (my$g=1; $g < @header; $g++){
			if(exists($genomes1{ $header[$g] } )){
				$genomes1{ $header[$g] } = 2;
				unless($q){print "\tList #1\t".$g."--".$header[$g]."\n"};
				push(@ind1,$g);
				push(@listnames,basename($list1));
			}elsif(exists($genomes2{ $header[$g] } )){
				$genomes2{ $header[$g] } = 2;
				unless($q){print "\tList #2\t".$g."--".$header[$g]."\n"};
				push(@ind2,$g);
				push(@listnames,basename($list2));
			}else{
				push(@others,$g);
				push(@listnames,'other');
			}
		}
		
		print OUT join("\t",@listnames)."\n"; #print row of list ids
		
		$c=1;
		unless($q){
			print "  Input matrix contains ". (@others + @ind1 + @ind2)." genomes.\n";
			print "  Found ".@ind1."/".(keys%genomes1)." genomes in list #1 and ".@ind2."/". (keys%genomes2) . " genomes in list #2.\n";
			if((@ind1+@ind2) < ((keys%genomes1)+(keys%genomes2))){
				print "  Missing genomes:\n";
				foreach my$mis (sort(keys%genomes1)){
					if($genomes1{$mis} == 1){ print "\t".$mis."\n"};
				}
				foreach my$mis2 (sort(keys%genomes2)){
					if($genomes2{$mis2} == 1){ print "\t".$mis2."\n"};
				}
			}
		}
	}else{

		my@c1=(0,0);
		my@c2=(0,0);
		chomp$i;
		my@adjs=split("\t",$i);
		
		## store matrix row ##
		$adjrows{$adjs[0]} = $i;
		
		## count orthologs ##
		foreach my$ii (@ind1){
			if($adjs[$ii] == 1 ){
				$c1[0]++;
			}else{$c1[1]++};
		}
		foreach my$io (@ind2){
			if($adjs[$io] == 1 ){
				$c2[0]++;
			}else{$c2[1]++};
		}
		
		## Run Fisher's Exact Test in R and get p-value ##
		$R->run( qq' fet = fisher.test(t(matrix(c($c1[0],$c1[1],$c2[0],$c2[1]),nrow=2))) ');
		my$fisher=$R->get('fet$p.value');
		
		## Store counts and p-value for each adjacency ##
		my$c1s = "(".join(",",@c1).")";
		my$c2s = "(".join(",",@c2).")";
		@{ $adjhash{ $j }} = ($adjs[0],$c1s,$c2s,$fisher);

		$j++;
	}

}
close IN;

## Benjamini-Hochberg P-value Adjustment ##
my@pvalues=();
foreach my$pp (natsort(keys%adjhash)){
	#print $pp."--";
	push(@pvalues, ${ $adjhash{ $pp} }[3]);
}

$R->set('ps', \@pvalues);
$R->run( qq' pa = p.adjust(ps, method="BH")' );
my$padj = $R->get('pa');

## Filter based on p-value and output list of passed to summary statistics file ## 
open STAT, ">$stats";
print STAT basename($in)."\t".basename($list1)."\t".basename($list2)."\tP.value\tBH-P.value\n";
my$pi=0;
for(my$pout=0; $pout < @pvalues; $pout++){
	if($pvalues[$pout] <= $p){
		print STAT @{$adjhash{$pout}}[0]."\t".@{$adjhash{$pout}}[1]."\t".@{$adjhash{$pout}}[2]."\t".sprintf("%.6f",@{$adjhash{$pout}}[3])."\t".sprintf("%.6f",$$padj[$pout])."\n";
		print OUT $adjrows{ @{$adjhash{$pout}}[0] } ."\n";
		
		$pi++;
	}
}
close STAT;
unless($q){print "  Found ".$pi." adjancies with p-values <= ".$p.".\n"};


##########################
sub HELP_MESSAGE { die "
.Description:	
   Takes binary table of specific 'adjaceny' presence/absence and compares presence between two lists of genomes by Fisher's Exact Test (using R).

.Usage: $0 -in in.txt -l1 list1.txt -l2 list2.txt -out out.txt
   
   [mandatory]
   -in 		<in.txt>	Input binary adjacency table, like the output from 'mauve.matrix2adj.pl'
   -l1		<list1.txt>	First list of genomes for comparison.
   -l2		<list2.txt>	Other list of genomes for comparison.
   -out		<out.txt>	New output table.
   -stats	<stats.txt>	New summary table of adjacencies with p-values below cutoff.
   
   [optional]
   -p		<0.05>		P-value cutoff (default 0.05).
   -c				Correct p-values for multiple testing.
   -q				Run quietly.

" }
