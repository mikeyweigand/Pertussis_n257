#!/bin/bash

if [[ "$1" == "" || "$1" == "-h" || "$1" == "-help" ]] ; then
   echo "
   This script will align a query genome with each reference genome in a list, individually, by many parallel runs of progressiveMauve.

   Usage: ./colinear-test-parallel.sh [query] [ref-list] [out-dir] [path-to-mauve]

   [mandatory]
   query	A closed genome sequence in Genbank(.gbk) or Fasta(.fasta) format.
   ref-list	List of paths to set of reference genomes in Genbank(.gbk) or Fasta(.fasta) format.
   out-dir	Path to the new output dir

   [optional]
   progressiveMauve   Full path to your progressiveMauve executable if not in \$PATH

   [dependencies]
   progressiveMauve
   parallel (https://www.gnu.org/software/parallel/)

   " >&2 ;
   exit 1 ;
fi ;

#where's your mauve?
if [[ "$4" == "" ]] ; then
  type progressiveMauve >/dev/null 2>&1 || { echo -e >&2 "\nERROR: Please indicate the full path to progressiveMauve. Run './colinear-test-parallel.sh -h' for details\n"; exit 1;}
else
  type $4 >/dev/null 2>&1 || { echo -e >&2 "\nERROR: '$4' not found. Please indicate the full path to progressiveMauve. Run './colinear-test-parallel.sh -h' for details\n"; exit 1;}
fi


## Define function for performing pairwise mauve alignments.
function collinear {

	name1=$(basename $1 | sed 's/.fasta//g' | sed 's/.gbk//g');
	name2=$(basename $2 | sed 's/.fasta//g' | sed 's/.gbk//g');
  #full1=$(find `pwd` -name "$(basename $1)");
  #full2=$(find `pwd` -name "$(basename $2)");
  full1=$(find `pwd` -path "*$1");
  full2=$(find `pwd` -path "*$2");


	echo $name1 vs $name2;
	d=$(date +%H%M%S%N);
	td=$3"/tmp_"$d;
	out=$3/$name1\-$name2.xmfa;
	log=$3/$name1\-$name2.log;

  echo -e "\t"$full1"\n\t"$full2;

	if [ ! -d "$3" ]; then
		mkdir $3;
	fi

	mkdir $td;
	cp $1 $td/;
	cp $2 $td/;

	date > $log;
	echo >> $log;
	$mymauve --output=$out --seed-weight=16 --hmm-identity=0.85 $td/$(basename $1) $td/$(basename $2) >> $log;

  #correct file paths in xmfa
  prev1="$td/$(basename $1)";
  prev2="$td/$(basename $2)";
  sed -i "s~$prev1~$full1~g" "$out";
  sed -i "s~$prev2~$full2~g" "$out";
  sed -i "s~$3~\.~g" "$out";

  gzip $out;
	echo >> $log;
	date >> $log;
	rm -r $td;
}
export -f collinear

## Run all pairwise alignments with parallel.
genome=$1;
reflist=$2;
outdir=$3;

if [ ! -d "$outdir" ]; then
	mkdir $outdir;
fi

mapfile -t refs < $reflist;
#parallel -j 12 -k echo $genome {} $outdir ::: ${refs[@]};
parallel -j 12 -k collinear $genome {} $outdir ::: ${refs[@]};
