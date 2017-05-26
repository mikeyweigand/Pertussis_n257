#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

&GetOptions(	'in=s' => \my$glist,		#
		'q' => \my$quiet,
		'h' => \my$h,
		'help' => \my$help,
		'tmp=s' => \my$tmp,
		'check' => \my$check,
		'dir=s' => \my$dir);		#
($glist and $tmp and $dir) or &HELP_MESSAGE;
if($help or $h){ &HELP_MESSAGE};

#read input list
my@genomes = ();
my$found=0;
open IN, "$glist";
while(my$i = <IN>){
	chomp$i;
	$i =~ s/^\.\///g; #remove preceeding "./" in file path, if present
	push(@genomes,$i);

	if($check){
		if( -e $i ){ $found++;
		}else{ print "\n\tERROR: Cannot find: ".$i."\n"};
	}

}
close IN;

#calculates total number of alignments
#my$numer = &factorial(scalar(@genomes));
#my$denom = 2*&factorial(scalar(@genomes)-2);
#my$numalign = $numer/$denom;
my$numalign = scalar(@genomes) * (scalar(@genomes)-1 ) / 2;

if($check){
	print "Found ".$found."/".scalar(@genomes)." files listed in '$glist'.\n";
	if( -d "$tmp"){print "WARNING: The directory '$tmp' already exists.\n"};
	if( -d "$dir"){print "WARNING: The directory '$dir' already exists.\n"};
	print "Total number of alignments = $numalign. No alignments were run.\n";
}

unless($check){
	#check for output directory
	unless( -d "$tmp") {
		system("mkdir $tmp");
	}
	unless($quiet){print "\nPerforming $numalign pairwise alignments...\n"};
	#prepare inputs and run pairwise-mauve-parallel.sh
	my@gen2 = @genomes;
	my$total = scalar(@genomes)-1;
	for(my$g=0; $g < scalar(@genomes)-1; $g++){
		shift@gen2;

		#save reference list
		my$tmpfile = basename($genomes[$g]);
		$tmpfile =~ s/\.\w+$/-list.txt/g;
		open TMP, ">$tmp/$tmpfile";
		print TMP join("\n",@gen2)."\n";
		close TMP;

		unless($quiet){ print "\n" . ($g+1) ."/".$total."\t". $genomes[$g]." vs ". $tmp."/".$tmpfile." (n=".scalar(@gen2).")\n"};

		#run parallel alignments
		system(" /home/yrh8/Documents/Pertussis_n257/Current/pairwise-mauve-parallel.sh $genomes[$g] $tmp/$tmpfile $dir");

	}
	unless($quiet){print "Done!\n\n"};
}
#######################
sub factorial {
	my$n = $_[0];
	my$fact = 1;
	for(my$o = 1; $o <= $n; $o++){
		$fact *= $o;
	}
	return( $fact )
}

sub HELP_MESSAGE { die "
.Description:
   Takes list of genome sequence file paths and computes all possible pairwise mauve alignments using 'pairwise-mauve-parallel.sh'.

.Usage: $0 -in -tmp -dir

   [mandatory]
	 -in	<in.txt>	List of paths to genome sequence files in genbank or fasta format.
	 -tmp	<path>		Directory for intermediate files.
	 -dir 	<path>		Directory for mauve alignment output files.

   [optional]
	 -q	Run quietly.
	 -check	Checks that all listed genome files and output directories exist without running alignments.

   [dependencies]
	 pairwise-mauve-parallel.sh.

" }
