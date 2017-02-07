#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;

&GetOptions(	'in=s' => \my$backbone,		#
		'gap=s' => \my$gap,		#
		'q' => \my$quiet,
		'out=s' => \my$out);		#
($backbone and $gap and $out) or &HELP_MESSAGE;

my($ext) = $backbone =~ /(\.xmfa[\.\w]+$)/;
my$inname= basename($backbone, $ext);# ".xmfa.backbone");
my@ids = split("-",$inname);

my$gapcount=0;
my$negcount=0;

open BB, "$backbone";
while(my$block = <BB>){
	unless($block =~ /seq/){
		chomp$block;
		my@sblock = split("\t",$block);
		$gapcount += &gapped(\@sblock, $gap);
		$negcount += &inversion(\@sblock);
	}
}
close BB;

open OUT, ">$out";
print OUT join("\t",(@ids,$gapcount,$negcount))."\n";
unless($quiet){ print join("\t",(@ids,$gapcount,$negcount))."\n"};
close OUT;

#####################
sub inversion {
	my@coords = @{$_[0]};
	my$neg=0;
	foreach my$c (@coords){
		if($c < 0){
			$neg=1;
			last;
		}
	}
	return( $neg );
}

sub gapped {
	my@coords = @{$_[0]};
	my$glen = $_[1];
	my@lens=();
	my$g=0;
	for(my$c=0; $c < scalar(@coords); $c+=2){
		my$len = abs($coords[$c] - $coords[$c+1]);
		push(@lens,$len);
	}
	my@sortlens = sort{$a <=> $b}@lens;
	if( ($sortlens[0] == 0) && ($sortlens[-1] >= $glen)){
		$g=1;
	}
	return( $g );
}

sub HELP_MESSAGE { die "
.Description:
   Checks a mauve alignment for gaps and inversions to determine if two genomes are collinear based on defined parameters.

.Usage: $0 -in -gap -out

   [mandatory]
	 -in	<in.backbone>	Backbone file output from progressiveMauve.
	 			(Assumes format: \"GenomeA-GenomeB.xmfa.backbone\")
	 -gap	<number>	Maximum allowable gap size in bp (eg. 1000).
	 -out	<out.txt>	Tab-separated output summary file.
	 			(Outputs: GenomeA GenomeB [gaps] [inversions])

   [optional]
	 -q		Run quietly.

   [dependencies]


" }
