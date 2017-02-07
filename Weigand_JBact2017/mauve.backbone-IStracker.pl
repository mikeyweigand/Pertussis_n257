#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Sort::Key::Natural qw(natsort);
use List::Util qw(max);

&GetOptions(	'bbone=s' => \my$bbone,		# progressiveMauve backbone file
		'xmfa=s' => \my$xmfa,		# progressiveMauve xmfa file
		'bls=s' => \my$blsd,		# directory of blast output files with IS-element coordinates
		'q' => \my$q,
		'h' => \my$h,
		'sum=s' => \my$summ,		# output optional summary table of IS-element counts in each genome
		'out=s' => \my$out,		# output matrix of matched IS-elements
		'blocks=s' => \my$bout);	# output block matrix

($bbone and $xmfa and $blsd and $out and $bout) or &HELP_MESSAGE;
if($h){ &HELP_MESSAGE};

## check that backbone and xmfa file are from the same progressiveMauve run ##
unless(basename($bbone, ".backbone") eq basename($xmfa)){
	print "\nERROR: These input files don't match:\n\t". basename($bbone)."\n\t".basename($xmfa)."\n\n";
	exit;
}

## get genome list and seq numbers from xmfa ##
my%xseqs=(); #hash of sequence number/index and file basename; ex: {'0'} = 'H866_CP011710'
my@genomes=();
my%pdlid=();
unless($q){ print "\n\tReading xmfa file..."};
open XMFA, "$xmfa";
while(my$x = <XMFA>){
	chomp$x;
	if($x =~ /^#Sequence\d+File/){
		my($xi) = $x =~ /(?<=Sequence)(\d+)/;	#xmfa 'sequence number'
		my($ext) = $x =~ /(\.\w+$)/;
		my@sx = split("\t",$x);
		$xseqs{ $xi } = basename($sx[1],$ext);	#sequence {number} = file basename
		push(@genomes,basename($sx[1],$ext));

		my($pdl) = basename($sx[1],$ext) =~ /([A-Z][0-9]{3})(=?[-_])/;
		#my($pdl) = basename($sx[1],$ext) =~ /(\w+)(=?[-_])/;
		$pdlid{ $pdl } = basename($sx[1],$ext);
		#print $pdl."\n";
	}
}

close XMFA;
unless($q){ print "done\n"};

#### store backbone file block coordinates ####
my%blocks=();
my@bbheader=("block");
my$bnum=sprintf("%06d",0);
unless($q){ print "\tReading backbone file..."};
open BBONE, "$bbone";
while(my$bb = <BBONE>){
	chomp$bb;
	my@sb=split("\t",$bb);
	#relabel backbone header
	if($bb =~ /^seq/){
		for (my$z=0; $z < @sb; $z+=2){
			my($bbhi) = $sb[$z] =~ /(?<=seq)(\d+)/;
			push(@bbheader, ($xseqs{$bbhi+1}."_leftend",$xseqs{$bbhi+1}."_rightend",$xseqs{$bbhi+1}."_length"));
		}
	}else{
		foreach my$g (natsort(keys%xseqs)){
			@{$blocks{$bnum}{ $xseqs{$g} }} = ($sb[$g*2-2],$sb[$g*2-1]); #{block number}{genome} = (leftend,rightend)
			#print "{".$bnum."}{".$xseqs{$g} ."} = (".$sb[$g*2-2].",".$sb[$g*2-1].")\n";
		}
		$bnum++;
	}
}
close BBONE;
unless($q){ print "done\n"};

#print new output matrix of block coordinates and lengths.
unless($q){ print "\tWriting new backbone file..."};

open BOUT, ">$bout";
print BOUT  join("\t",@bbheader)."\n";
foreach my$bo1 (natsort(keys%blocks)){
	print BOUT "block_".sprintf("%06d",$bo1);
	foreach my$bo2 (natsort(keys%xseqs)){
		print BOUT "\t".join("\t",@{$blocks{$bo1}{$xseqs{$bo2}} });
		print BOUT "\t".abs( $blocks{$bo1}{$xseqs{$bo2}}[1] - $blocks{$bo1}{$xseqs{$bo2}}[0]);

		#save all block lengths.
		push( @{$blocks{$bo1}{lengths}}, abs($blocks{$bo1}{$xseqs{$bo2}}[1] - $blocks{$bo1}{$xseqs{$bo2}}[0]) );

	}
	print BOUT "\n";
}
close BOUT;

unless($q){ print "done\n"};


#### IS-element blast results ####
my%iselements=(); #stores all elements per genome
unless($q){ print "\tReading IS-element insertions..."};
while(my$bls = <"$blsd/*">){
	open BLAST, "$bls";
	#my($id) = $bls =~ /([A-Z][0-9]{3})/;
	my($id) = $bls =~ /(\w+)(=?-IS)/;
	#print $id . "\n";

	#confirm that genome is in alignment before inporting blast hits (20160419) ***only works for PDL strain names***
	if(exists($pdlid{ $id })){
	#	print $id."\n";


	while(my$hit = <BLAST>){
		chomp$hit;
		my@shit = split("\t",$hit);
		my$coords=$shit[-2]."-".$shit[-1];

		##establish object structure##
		unless(exists( $iselements{ $id }{ $shit[0]} )){
			@{ $iselements{ $id }{ $shit[0] }{gcoords} } = (); #genome coordinate array
			@{ $iselements{ $id }{ $shit[0] }{blocks} } = (); #matched block number array
			@{ $iselements{ $id }{ $shit[0] }{bcoord} } = (); #relative block coordinate array
			@{ $iselements{ $id }{ $shit[0] }{dir} } = (); #relative direction array
			@{ $iselements{ $id }{ $shit[0] }{dup} } = (); #duplication status array
			$iselements{ $id }{ $shit[0] }{dupcount}  = 0; #count of detected duplications
		}

		#check for duplicate/neighboring insertion#
		my$dup="unknown";
		my$dcount=0;
		foreach my$dd (@{ $iselements{ $id }{ $shit[0] }{gcoords} }){
			my@cc = split("-",$dd);
			if(abs($shit[-2] - $cc[1]) <= 10 || abs($shit[-1] - $cc[0]) <= 10){ #within 10bp of previous match
				$dup="yes";
				@{ $iselements{ $id }{ $shit[0] }{dup} }[$dcount] = "yes"; #changes duplication status of existing/matched entry
				$iselements{ $id }{ $shit[0] }{ dupcount }++;
			}
			$dcount++;
		}

		push( @{ $iselements{ $id }{ $shit[0]}{gcoords} },$coords); #stores genome coordinates
		push( @{ $iselements{ $id }{ $shit[0] }{dup} },$dup);
		#print "\t".$shit[0]."\t".$coords."\t".$dup."\t||\t";

		#identify block number containing IS-element
		my$mid = (abs($shit[-1]) + abs($shit[-2]))/2;
		my$nomatch=0;
		foreach my$bn (keys%blocks){
			my@bc = @{ $blocks{ $bn }{ $pdlid{$id} } };
			if( abs($bc[0]) <= $mid && $mid <= abs($bc[1]) ){

				#store matched block number
				push( @{ $iselements{ $id }{ $shit[0]}{blocks} }, sprintf("%06d",$bn)); #stores matched block number

				#store matched block relative coordinate and orientation
				my@rpos = REL_COORD($shit[-2],$shit[-1],$bc[0],$bc[1]);
				push( @{ $iselements{ $id }{ $shit[0]}{bcoord} },$rpos[0]); #stores relative position
				push( @{ $iselements{ $id }{ $shit[0]}{dir} },$rpos[1]); #stores relative direction

				#print sprintf("%06d",$bn)."\t".$rpos[0]."\t".$rpos[1];
				last;
			}else{
				$nomatch++;
			}
		}
		if( $nomatch == scalar(keys%blocks)){
			push( @{ $iselements{ $id }{ $shit[0]}{blocks} }, "NONE"); #stores matched block number
			push( @{ $iselements{ $id }{ $shit[0]}{bcoord} },"0"); #no block match
			push( @{ $iselements{ $id }{ $shit[0]}{dir} },"NA"); #no block match
			#print "NONE\tNA\tNA";
		}
		#print "\n";
	}
	}
}
unless($q){ print "done\n"};


### create non-redundant list based on block number, orientation, and block-relative position ###
#{IS481}{block#}{FWD/REV}{rel-position}{dup}
my%islist=(); #stores non-redundant list of 'orthologous/matched' elements
my%isout=(); #stores matched elements for simpler output format
my$isid='';
my$wtf=0;
unless($q){ print "\tMatching IS-element insertions..."};
foreach my$gen (natsort(keys%iselements)){
	#print $gen."\n"; #. scalar(keys$iselements{$gen})."\n";

	foreach my$ise (natsort(keys$iselements{$gen})){
		#print "\t".$ise."\t". scalar(@{ $iselements{ $gen }{ $ise }{gcoords} })."\n";

		my@ag = @{ $iselements{ $gen }{ $ise }{gcoords} }; #genome coordinate array
		my@abb = @{ $iselements{ $gen }{ $ise }{blocks} }; #matched block number array
		my@abc = @{ $iselements{ $gen }{ $ise }{bcoord} }; #relative block coordinate array
		my@ard = @{ $iselements{ $gen }{ $ise }{dir} }; #relative direction array
		my@adp = @{ $iselements{ $gen }{ $ise }{dup} }; #duplication status array

		for (my$ii=0; $ii < scalar(@ag); $ii++){
			#check for perfect match
			if( exists( $islist{$ise}{$abb[$ii]}{$ard[$ii]}{$abc[$ii]} )){

				push( @{ $islist{$ise}{$abb[$ii]}{$ard[$ii]}{$abc[$ii]} }, [$pdlid{ $gen },$ag[$ii]] );

				$isid=$abb[$ii]."_".$abc[$ii]."_".$ise."_".$ard[$ii];
				$isout{ $isid }{$pdlid{ $gen }} = $ag[$ii];

				if($adp[$ii] =~ /[Yy]/){
					$isout{$isid}{dup} = $adp[$ii]; #over-write dup status if 'yes'
				}

			}else{
				#check for imperfect match (+/-20bp)
				my$counter=0;
				foreach my$check ( keys$islist{$ise}{$abb[$ii]}{$ard[$ii]} ){

					if( abs($check - $abc[$ii]) <= 20){
						push( @{ $islist{$ise}{$abb[$ii]}{$ard[$ii]}{$check} }, [$pdlid{ $gen },$ag[$ii]]);

						$isid=$abb[$ii]."_".$check."_".$ise."_".$ard[$ii];
						$isout{ $isid }{$pdlid{ $gen }} = $ag[$ii];

						if($adp[$ii] =~ /[Yy]/){
							$isout{$isid}{dup} = $adp[$ii]; #over-write dup status if 'yes'
						}
						last;

					}else{ $counter++};
				}

				#store 'new' match if none found
				if( $counter == keys$islist{$ise}{$abb[$ii]}{$ard[$ii]}){
					@{ $islist{$ise}{$abb[$ii]}{$ard[$ii]}{$abc[$ii]} } = ();
					push( @{ $islist{$ise}{$abb[$ii]}{$ard[$ii]}{$abc[$ii]} }, [$pdlid{ $gen },$ag[$ii]] );

					$isid=$abb[$ii]."_".$abc[$ii]."_".$ise."_".$ard[$ii];
					$isout{ $isid }{$pdlid{ $gen }} = $ag[$ii];
					$isout{ $isid }{dup} = $adp[$ii]; #duplication status
					$wtf++;
				}
			}
		}

	}

}
unless($q){ print "done\n"};

#print summary of matched IS-elements in each genome to STOUT
unless($q){print "\n\tThere are " .scalar(@genomes)." genomes, with ".scalar(keys%isout)." unique IS-element insertions:\n\n"};
my%isout2=();

foreach my$o (natsort(keys%iselements)){
	unless($q){ print "\t". $o.",".$pdlid{$o}."\n"};
	foreach my$k (natsort(keys$iselements{$o})){
		unless($q){ print "\t\t\t".$k."\t". scalar(@{ $iselements{$o}{$k}{gcoords}})."\t(". $iselements{$o}{$k}{dupcount} ." duplications)\n"};
		$isout2{$k}{$pdlid{$o}} = scalar(@{ $iselements{$o}{$k}{gcoords}} );
	}
	unless($q){ print "\n" };
}

#print new output matrix of matched IS-elements.
unless($q){ print "\tWriting new matrix of matched IS-element insertions..."};
open OUT, ">$out";
print OUT "IS-element\tMax-block-len\tDuplication\tCount\t". join("\t",@genomes)."\n";
foreach my$is1 (natsort(keys%isout)){
	print OUT $is1;
	#print $is1."\n";

	if($is1 =~ /^NONE/){
		print OUT "\t0";
	}else{
		my($bo3) = $is1 =~ /^([0-9]{6})(?=_)/; #block number
		print OUT "\t".max(@{$blocks{ $bo3}{lengths}}); #maximum block length
	}
	print OUT "\t".$isout{$is1}{dup};

	my$ccount = scalar(keys%{$isout{$is1}}) - 1;
	print OUT "\t". $ccount;

	foreach my$is2 (@genomes){
		if(exists $isout{$is1}{$is2}){ print OUT "\t".$isout{$is1}{$is2};
		}else{ print OUT "\tNA";
		}

	}
	print OUT "\n";
}

unless($q){ print "done\n"};
#print optional table of IS-element counts in each genome
if($summ){
	open SUM, ">$summ";
	print SUM "IS-element\t". join("\t",natsort(@genomes))."\n";
	foreach my$is3 (natsort(keys%isout2)){
		print SUM $is3;
		foreach my$go1 (natsort(@genomes)){
			print SUM "\t". $isout2{$is3}{$go1};
		}
		print SUM "\n";
	}
}

##########################
sub REL_COORD {
	my($is1,$is2,$mbc1,$mbc2) = @_;
	my@rel=();

	unless($mbc1 < 0){
		if($is1 < $is2){
			@rel=(($is1-$mbc1+1),'FWD');
		}else{
			@rel=(($is2-$mbc1+1),'REV');
		}
	}else{
		if($is1 < $is2){
			@rel=((abs($mbc2)-$is2+1),'REV');
		}else{
			@rel=((abs($mbc2)-$is1+1),'FWD');
		}
	}
	return @rel;
}

sub HELP_MESSAGE { die "
.Description:
   Matches IS-element insertions across genomes using progressiveMauve (xmfa and backbone) and BLASTn outputs.

.Usage: $0 -bbone in.bbone -xmfa in.xmfa -bls blastout/ -out matrix.txt -blocks blocks.txt -sum sum.txt

   [mandatory]
   -bbone	<in.bbone>	Input backbone file from progressiveMauve (must match xmfa file).
   -xmfa	<in.xmfa>	Input xmfa file from progressiveMauve (must match backbone file).
   -bls		<dir>		Directory of tabular BLASTn output files.
   				Expected outfmt: 'qseqid sseqid qlen length qstart qend sstart send'
   -out		<matrix.txt>	New matrix of matched IS-element insertions.
   -blocks	<blocks.txt>	New backbone table of block coordinates and lengths.

   [optional]
   -q				Run quietly.
   -h				This helpful message.
   -sum		<sum.txt>	Optional summary table of total IS-element counts per genome.

 WARNING: Written to specifically recognize PDL strain-ids and GenBank accession numbers, eg. A001_CP000001. 
" }
