#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

#description: create signed permutation matrix from progressiveMauve *.backbone file.

&GetOptions(	
		'nr=s' => \my$nr,		# permutation matrix of non-redundant subset
		'list=s' => \my$list,		# full list of all genomes
		'out=s' => \my$out);		# new permutation matrix with all genomes

($nr and $list and $out) or &HELP_MESSAGE;

#store genome id and group name
my%hash=();
open LIST, "$list";
while(my$i = <LIST>){
	#print $i;
	chomp$i;
	my$name = basename($i);
	$name =~ s/\.[a-z]+$//;
	#print "\t".$name."\t";
	
	my($group) = $i =~ /(?<=\/)([\w-]+)(?=\/)/;
	#print $group ."\n";
	
	if($group =~ /orphan|vaccine|ptxP/){
		#print "\t".$name."\t".$name."\n";
		$hash{ $name } = $name;
	}else{
		$hash{ $name } = $group;
		#print $group. "::".$name."\n";
	}
}
close LIST;

#store set of non-redundant, reference genomes
my%ref=();
open NR, "$nr";
my@nr_matrix=();
while(my$n = <NR>){
	chomp$n;
	push(@nr_matrix,$n);
}
close NR;

for( my$c=0; $c < @nr_matrix; $c+=2){
	$nr_matrix[$c] =~ s/^>//;
	#print $nr_matrix[$c] . "\t". $hash{$nr_matrix[$c]} ."\n";
	$ref{ $hash{$nr_matrix[$c]} } = $nr_matrix[$c+1];
}

#apply reference permutation matrix to all genomes
open OUT, ">$out";
my$m=0;
foreach my$go (sort(keys%hash)){
	if(exists($ref{ $hash{$go} })){
		print OUT ">".$go."\n".$ref{ $hash{$go} }."\n";
		$m++;
	}else{
		print "\tERROR: $go HAS NO MATCH!\n"
	}
}
close OUT;

print "\n\t$m / ".scalar(keys%hash)." genomes matched to ". @nr_matrix/2 ." non-redundant references.\n\n";

sub HELP_MESSAGE { die "
.Description:
   Transfers signed permutation matrix across co-linear genomes, for input to MLGO.
   
.Usage: $0 -nr <matrix.txt> -list <list.txt> -out <out.txt>
   
   [mandatory options]
   -nr		Input permutation matrix of Non-redundant subset.
   -list	All genomes, in group sub-directories.
   		Ex: /CDC010/H806_CP000000.gbk
   -out		Filename of new output matrix.

" }
