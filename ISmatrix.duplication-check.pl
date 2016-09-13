#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Sort::Key::Natural qw(natsort);
use List::Util qw(max);

#default paramters:
my$flank = 10;

&GetOptions(	'bls=s' => \my$blsd,		# directory of blast output files with IS-element coordinates
		'q' => \my$q,
		'h' => \my$h,
		'flank=s' => $flank,
		'out=s' => \my$out,		# new output matrix of IS-elements
		'ogs=s' => \my$matrix);	# existing matrix of IS-elements

($matrix and $blsd and $out) or &HELP_MESSAGE;
if($h){ &HELP_MESSAGE};

#### Grab genome list from matrix header ####
my$head = qx(head -1 $matrix);
chomp$head;
my@header = split(/\s+/,$head);
my%pdlid = ();
for(my$n = 3;$n < scalar@header; $n++){
	my($pid) = $header[$n] =~ /^([A-Z][0-9]{3})/;
	my($genbank) = $header[$n] =~ /(CP[0-9]{6})$/;
	$pdlid{$pid}=$genbank;
#	print $pid."--".$genbank."\n";
}

#### IS-element blast results ####
my%iselements=(); #stores all elements per genome
my%dup=();
my@dall=();
#unless($q){ print "\tReading IS-element insertions...\n"};
while(my$bls = <"$blsd/*">){
	open BLAST, "$bls";
	#my($id) = $bls =~ /([A-Z][0-9]{3})/;
	my($id) = $bls =~ /(\w+)(=?-IS)/;
	#print $id . "\n";

	#confirm that genome is in alignment before inporting blast hits (20160419) ***only works for PDL strain names***
	if(exists($pdlid{ $id })){

		while(my$hit = <BLAST>){
			chomp$hit;
			my@shit = split("\t",$hit);
			my$coords=$shit[-2]."-".$shit[-1];
			$shit[0] =~ s/v\d//g;
			##establish object structure##
			unless(exists( $iselements{ $id }{ $shit[0]} )){
				#print $id."--".$shit[0]."\n";
				@{ $iselements{ $id }{ $shit[0] }{gcoords} } = (); #genome coordinate array
				#@{ $iselements{ $id }{ $shit[0] }{blocks} } = (); #matched block number array
				#@{ $iselements{ $id }{ $shit[0] }{bcoord} } = (); #relative block coordinate array
				#@{ $iselements{ $id }{ $shit[0] }{dir} } = (); #relative direction array
				#@{ $iselements{ $id }{ $shit[0] }{dup} } = (); #duplication status array
				#$iselements{ $id }{ $shit[0] }{dupcount}  = 0; #count of detected duplications
			}

			#check for duplicate/neighboring insertion#
			my$dup="unknown";
			my$dcount=0;
			foreach my$dd (@{ $iselements{ $id }{ $shit[0] }{gcoords} }){
				my@cc = split("-",$dd);
				if(abs($shit[-2] - $cc[1]) <= $flank || abs($shit[-1] - $cc[0]) <= $flank){ #within 10bp of previous match
					#$dup="yes";
					#@{ $iselements{ $id }{ $shit[0] }{dup} }[$dcount] = "yes"; #changes duplication status of existing/matched entry
					#$iselements{ $id }{ $shit[0] }{ dupcount }++;

					my$cd = $shit[-2]."-".$shit[-1];

					if(exists( $dup{ $id }{ $shit[0] }{ $cd })){
						push( @{ $dup{ $id }{ $shit[0] }{ $cd } },$dd);
					}elsif(exists( $dup{ $id }{ $shit[0] }{ $dd })){
						push( @{ $dup{ $id }{ $shit[0] }{ $dd } },$cd);
					}else{
						@{ $dup{ $id }{ $shit[0] }{ $dd } } = ($cd);
					}
					#$dup{ $id }{ $shit[0] }{ $cd } = $dd;
					#print $id."--".$shit[0].":".$cd."-". $dd."\n";
				}
				$dcount++;
			}
			#print $id."--".$shit[0]."--".$coords."\n";
			push( @{ $iselements{ $id }{ $shit[0]}{gcoords} },$coords); #stores genome coordinates
			#push( @{ $iselements{ $id }{ $shit[0] }{dup} },$dup);
			#print "\t".$shit[0]."\t".$coords."\t".$dup."\t||\t";
		}
		#print $id."\n";
		foreach my$o1 (keys$dup{ $id }){
				#print $o1."\n";
			foreach my$o2 (natsort(keys$dup{$id}{$o1})){
				my@o3 = @{ $dup{$id}{$o1}{$o2} };
				push(@o3,$o2);
				#print "\t\t".join("||",sort@o3)."\t\n";

				foreach my$dl (@o3){
					my$dsave = $id."_".$pdlid{$id}.".".$o1.".".$dl;
					push(@dall,$dsave);
					#print"\t". $dsave."\n";
				}
				for (my$dlc=0;$dlc<scalar(@o3);$dlc++){
					my@dlout = @o3;
					splice(@dlout,$dlc,1);
					#print "\t\t\t\t\t\t\t\t". $o3[$dlc]."\t\t".join("--",@dlout)."\n";
					@{ $dup{$id}{$o1}{$o3[$dlc]} } = @dlout;

				}

			}

		}


	}
	close BLAST;
}


### open current IS matrix ###
open MX, "$matrix";
open OUT, ">$out";
my%fixed=();

while(my$m = <MX>){

	if($m =~ /^tag/){
		print OUT $m;
	}else{
		my$dcount=0;
		foreach my$dcheck (@dall){	#check for presence of duplicated IS-insertion
			if($m =~ /$dcheck/){
				last;
			}else{
				$dcount++;
			}
		}
		#print scalar(@dall)."--". $dcount."\n";;
		if($dcount == scalar(@dall)){
			print OUT $m;	#print non-duplicated sites
		}else{
			chomp$m;
			my@sm = split("\t",$m);
			#print $sm[0]."\t".$sm[1]."\n";

			for(my$mt=3;$mt < scalar(@sm); $mt++){
				if(length($sm[$mt]) > 1){ #make sure value isn't "-"
					#print $sm[$mt]."\t\t";
					my@ssm = split(",",$sm[$mt]);
					my@sssm = split(/\./,$ssm[0]);
					my$isid= join("\.",($sssm[0],$sssm[1])).".";
					my$mid = substr($sm[$mt],0,4); #."\t";
					#print $mid."\t". $sssm[-1]."\t";

					if( exists($dup{ $mid }{$sssm[1]}{$sssm[-1]}) ){
						#print "yes\t".$isid.$sssm[-1];
						my$replace=$isid.$sssm[-1];
						foreach my$kout (@{ $dup{ $mid }{$sssm[1]}{$sssm[-1]} }){
							$replace.= ",". $isid.$kout;
						}
						$sm[2] = "yes";
						#print"--". $replace."\n";
						$sm[$mt] = $replace;
					}else{ }
				}
			}

			#check if any sites are repeated before printing to output
			my$prev=0;
			for(my$ch=3;$ch < scalar(@sm);$ch++){
				my@sch=split(",",$sm[$ch]);
				foreach my$ch2 (@sch){
					if(exists($fixed{$ch2})){
						$prev=1;
					}else{
						$fixed{$ch2}=1;
					}
				}
			}
			if($prev == 0){
				print OUT join("\t",@sm)."\n";
			}
		}
	}
}


sub HELP_MESSAGE { die "
.Description:
   Performs one final check of IS element duplications.

.Usage: $0 -ogs -bls -out

   [mandatory]
   -ogs	<in.ogs>	Matrix of IS matches, from 'ogs.rb' or 'ISmatrix.add-strain-spec.pl'.
   -bls	<dir>		Directory of tabular BLASTn output files with IS hits.
   			Expected outfmt: '6 qseqid sseqid qlen length qstart qend sstart send'
   -out	<out.ogs>	New matrix of IS matches.

   [optional]
	-flank	<number>	Flanking distance (bp) for considering duplication. (Default = 10)
	-h			This helpful message.

WARNING: Written to specifically recognize name formats containing PDL strain-ids and GenBank accession numbers, eg. 'A001_CP000001'.

" }
