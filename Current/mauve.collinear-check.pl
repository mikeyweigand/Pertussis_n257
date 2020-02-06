#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use POSIX qw(ceil);

&GetOptions(	'in=s' => \my$backbone,		#
		'gap=s' => \my$gap,		#
		'q' => \my$quiet,
		'coords=s' => \my$coordtable,
		'out=s' => \my$out);		#
($backbone and $gap and $out) or &HELP_MESSAGE;

my($ext) = $backbone =~ /(\.xmfa[\.\w]+$)/;
my$inname= basename($backbone, $ext);# ".xmfa.backbone");
my@ids = split("-",$inname);

my$gapcount=0;
my$invertcount=0;
my%coords=();
my%gaps=();
my$ind=0;

if($coordtable){ open CT, ">$coordtable"};

open BB, "$backbone";
while(my$block = <BB>){
	if($block =~ /seq/){
		if($coordtable){ print CT $block };
	}else{
		chomp$block;
		my@sblock = split("\t",$block);
		my$isgap = &gapped(\@sblock, $gap);

		if($isgap == 1){
			$gapcount += $isgap;
			$gaps{ $gapcount } = \@sblock;
		}elsif($isgap == 0){
			$coords{ abs($sblock[0]) }[0] = $ind;
			$coords{ abs($sblock[0]) }[1] = $sblock[0]; #stores coordinate
			$coords{ abs($sblock[0]) }[3] = \@sblock;
			if($sblock[0] > 0){
				$coords{ abs($sblock[0]) }[2] = 0; #positive, sign
			}else{
				$coords{ abs($sblock[0]) }[2] = 1;	#negative
			}
			$ind++;
		}
	}
}
close BB;

my@sorted = sort{$a <=> $b}(keys(%coords)); #sorts by abs value of coord
my@outcoords = @{ $coords{ $sorted[0] }[3] };
my@gblock=();
for(my$s=1	; $s < @sorted; $s++){
	if( $coords{$sorted[$s]}[2] == $coords{$sorted[$s-1]}[2] ){	#same sign?
		if ( abs($coords{$sorted[$s]}[0] - $coords{$sorted[$s-1]}[0]) == 1 ){	#indices in order?

			@gblock = &gapcheck( \@outcoords, \@{$coords{$sorted[$s]}[3]}, \%gaps ); #check if disrupted by gap

			if(@gblock >= 1){
				if($coordtable){ print CT join("\t",@outcoords)."\n" };
				@outcoords = @{$coords{$sorted[$s]}[3]};

			}else{
				#combine with preceeding block
				$outcoords[1] = @{$coords{$sorted[$s]}[3]}[1];
				if( $coords{$sorted[$s]}[2] == 0 ){ #positive
					$outcoords[3] = @{$coords{$sorted[$s]}[3]}[3];
				}else{ #negative
					$outcoords[2] = @{$coords{$sorted[$s]}[3]}[2];
				}
			}
		}else{
			#rearrangement, start new block
			if($coordtable){ print CT join("\t",@outcoords)."\n" };
			@outcoords = @{$coords{$sorted[$s]}[3]};
			$invertcount++;
		}
	}else{
		#inversion, start new block
		if($coordtable){ print CT join("\t",@outcoords)."\n" };
		@outcoords = @{$coords{$sorted[$s]}[3]};
		$invertcount++;
	}
}
if($coordtable){
	print CT join("\t",@outcoords)."\n";
	foreach my$gg (sort{$a <=> $b}(keys(%gaps))){
		print CT join("\t",@{$gaps{$gg}})."\n";
	}
}

my$negcount = ceil($invertcount / 2);
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
	if(($sortlens[0] == 0) && ($sortlens[-1] >= $glen)){
		$g=1;
	}elsif(($sortlens[0] == 0) && ($sortlens[-1] < $glen)){
		$g=2;
	}
	return( $g );
}

sub gapcheck {
	my@prev = @{$_[0]};
	my@next = @{$_[1]};
	my%gaps = %{$_[2]};
	my@out = ();
	foreach my$g (sort{$a <=> $b}(keys(%gaps))){
		my@gblock = @{$gaps{$g}};
		if(($gblock[0] > $prev[1]) && ($gblock[1] < $next[0])){
			@out = @gblock;
		}elsif(($gblock[2] > $prev[3]) && ($gblock[3] < $next[2])){
			@out = @gblock;
		}
	}
	return( @out );

}

sub HELP_MESSAGE { die "
.Description:
   Checks a mauve alignment for gaps and inversions to determine if two genomes are collinear based on defined parameters.

.Usage: $0 -in -gap -out

   [mandatory]
	 -in	<in.backbone>	Backbone file output from progressiveMauve.
	 			(Assumes format: \"GenomeA-GenomeB.xmfa.backbone[.clean]\" with no other \"-\" in the filename)
	 -gap	<number>	Maximum allowable gap size in bp (eg. 1000).
	 -out	<out.txt>	Tab-separated output summary file.
	 			(Outputs: GenomeA GenomeB [gaps] [inversions])

   [optional]
	 -q				Run quietly.
	 -coords	<coords.txt>	Simplified backbone file of gap and inversion coordinates.

   [dependencies]


" }
