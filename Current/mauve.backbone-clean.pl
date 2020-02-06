#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Bio::SeqIO;

&GetOptions(	'in=s' => \my$backbone,		#
		'query=s' => \my$query,		#
		'q' => \my$quiet,
		'h' => \my$h,
		'genomes=s' => \my$genomes,
		'out=s' => \my$out);		#
($backbone and $query and $genomes and $out) or &HELP_MESSAGE;
if($h){ &HELP_MESSAGE };

my$inname= basename($backbone, ".xmfa.backbone");
my@ids = split("-",$inname);
my$genomeA = &getgrep( $ids[0], $genomes);
my$genomeB = &getgrep( $ids[1], $genomes);

my$clean=0;

open BB, "$backbone";
open OUT, ">$out";
while(my$block = <BB>){
	if($block =~ /seq/){
		print OUT $block;
	}else{
		chomp$block;
		my@sblock = split("\t",$block);
		if( &gapcheck(\@sblock,$query,$genomeA,$genomeB) == 1){
			print OUT $block ."\n";
		}else{ $clean++ };
	}
}
close BB;
close OUT;

unless($quiet){
	print "\nRemoved ".$clean." gaps from '".basename($backbone)."' based on matches to '".basename($query)."'.\n\n";
}
#clean-up any temporary fasta files created for blasting
foreach my$g ($genomeA,$genomeB){
	if($g =~ /tmp.fasta$/){
		system( "rm $g" );
	}
}

##############################
sub gapcheck {
	my@coords = @{$_[0]};
	my%coordhash = map { $_ => 1 } @coords;
	if(exists($coordhash{ 0 })){
		my$r=0;
		if($coords[0] == 0){ #blast genomeB
			$r = &doblast($_[1],$_[3],$coords[2],$coords[3]);
		}elsif($coords[2] == 0){ #blast genomeA
			$r = &doblast($_[1],$_[2],$coords[0],$coords[1]);
		}
		return($r);

	}else{ return(1) };
}

sub doblast {
	my($q, $sub, $c1, $c2) = @_;
	my$blast = qx( blastn -query $q -subject $sub -subject_loc $c1\-$c2 -outfmt \'6 sstart send\' -qcov_hsp_perc 50 -perc_identity 90 );
	chomp$blast;
	$blast =~ s/\n/\t/g;
	my@hits = sort{$a <=> $b}(split("\t",$blast));
	if( scalar(@hits) < 2){ #make sure there are hits
		return(1);
	}elsif((($hits[0] - $c1) >= 100) || (($c2 - $hits[-1]) >= 100)){ #matches cannot be >100bp inside gap
		return(1);
	}elsif( scalar(@hits) > 2){
		my$r=0;
		for(my$h = 2; $h < scalar(@hits); $h+=2){
			if(( $hits[$h] - $hits[$h-1] ) >= 50 ){ #max allowable 50bp between
				$r=1;
				last;
			}
		}
		return($r);
	}else{ return(0) };
}

sub getgrep {
	my($pdl) = $_[0] =~ /([A-Z][0-9]{3})/;
	my$match = qx( grep $pdl $_[1] );
	chomp$match;
	return( &checkfasta($match) );
}

sub checkfasta {
	if($_[0] =~ /fasta/){
		return ($_[0]);
	}else{ #convert genbank to temp fasta for blasting
		my$gbk = Bio::SeqIO->new(-file => "$_[0]", -format => 'genbank');
		my$tempfasta = basename($_[0],".gbk"). ".tmp.fasta";
		my$out = Bio::SeqIO->new(-file => ">$tempfasta", -format => 'fasta');
		while(my$seq = $gbk->next_seq){
			$out->write_seq($seq);
		}
		return($tempfasta);
	}
}


sub HELP_MESSAGE { die "
.Description:
   Checks a pairwise mauve alignment backbone file and removes gaps that correspond to BLASTn alignment with a given query (eg. IS481). Removes gaps if [1] query alignment coordinates are < 100bp within gap coordinates, and [2] multiple query matches are spaced < 50bp apart, where applicable.

.Usage: $0 -in -gap -out

   [mandatory]
	 -in		<in.backbone>	Backbone file output from progressiveMauve.
	 				(Assumes format: \"GenomeA-GenomeB.xmfa.backbone\")
	 -query		<query.fasta>	Fasta sequence for BLASTn.
	 -genomes 	<list.txt>	List of file paths to genome sequences in alignment.
	 				(Assumes basename consistency: GenomeA.gbk, GenomeB.fasta)
	 -out		<out.txt>	New backbone file with matched gaps removed.


   [optional]
	 -q		Run quietly.
	 -h		This help message.

   [dependencies]
	 BioPerl	(Bio::SeqIO)
	 blastn		(Must be in your \$PATH)


" }
