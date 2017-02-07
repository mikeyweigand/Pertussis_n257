#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Statistics::R;

#Default MCL parameters
my$inflation=2;

&GetOptions(	'in=s' => \my$check,		#
		'q' => \my$quiet,
		'S' => \my$singletons,
		'I=s' => \$inflation,
		'plot' => \my$plot,
		'out=s' => \my$out);		#
($check and $out) or &HELP_MESSAGE;

my%genomes=();
open CHECK, "$check";
open TMP, ">mcl-input.tmp";
while(my$c = <CHECK>){
	chomp$c;
	my@sc = split("\t",$c);
	$genomes{$sc[0]}=1;
	$genomes{$sc[1]}=1;
	if( ($sc[2] == 0) && ($sc[3] == 0)){ #no gaps, no inversions
		print TMP $sc[0]."\t".$sc[1]."\n";
	}
}
close CHECK;
close TMP;
#cluster matches using mcl
if($quiet){
	system( "mcl mcl-input.tmp -V all --abc -I $inflation -o mcl-output.tmp");
}else{
	system( "mcl mcl-input.tmp --abc -I $inflation -o mcl-output.tmp");
}

#add cluster names and counts
open MOUT, "mcl-output.tmp";
open COUT, ">$out";
#my$grouplist = basename($out, "txt")."groups.txt";
my$grouplist = $out;
$grouplist =~ s/txt$/groups.txt/;
open GOUT, ">$grouplist";
my$cnum=1;
my%hash = ();
while(my$m = <MOUT>){
	chomp$m;
	my@sm = split("\t",$m);
	my$cname = "Cluster-".$cnum;
	print COUT $cname."\t".scalar@sm."\t".$m."\n";
	#save cluster count for barplot
	$hash{$cname} = scalar@sm;
	#output list of cluster name and genome
	foreach my$glist (@sm){ print GOUT $cname."\t".$glist."\n"};

	$cnum++;
}
close MOUT;
close COUT;
system( "rm mcl-input.tmp mcl-output.tmp" );

#append unclustered genomes ('singletons')
if($singletons){
	my$scount=0;
	open OUT, ">>$out";
	foreach my$g (sort(keys%genomes)){
		if(single($g, $out) == 0){
			print OUT "Singleton\t1\t". $g ."\n";
			print GOUT "Singleton\t". $g."\n";
			$scount++;
		}
	}
	$hash{"Singletons"} = $scount;
	unless($quiet){ print $scount ." singletons found.\n\n" };
	close OUT;
	close GOUT;
}

#draw optional barplot
if($plot){
	my@k = sort{$hash{$b} <=> $hash{$a} }keys(%hash);
	my@v = @hash{@k};
	#my$plotout = basename($out, "txt")."png";
	my$plotout = $out;
	$plotout =~ s/txt$/png/;
	my$R = Statistics::R->new();
	$R->set('counts', \@v );
	$R->set('groups', \@k );
	$R->run( qq' png("$plotout",width=800,height=480) ');
	$R->run( qq' barplot(counts, names.arg = groups,ylab="Number of genomes",main = "$out",las=2) ');
	$R->run( qq' dev.off() ' );
}

########################
sub single {
	my$genome = $_[0];
	my$mcl = $_[1];
	my$match = qx( grep -c $genome $mcl );
	chomp$match;
	return( $match );
}

sub HELP_MESSAGE { die "
.Description:
   Takes concatenated output from 'mauve.collinear-check.pl' and groups collinear genomes using Markov Clustering Algorithm (MCL).

.Usage: $0 -in -out

   [mandatory]
	 -in	<in.txt>	Concatenated output from 'mauve.collinear-check.pl'.
	 			(Tab-separated format: GenomeA GenomeB [gaps] [inversions])
	 -out	<out.txt>	Output file of clusters from MCL.

   [optional]
	 -q	Run quietly.
	 -I	Set MCL 'inflation' parameter (default = 2.0).
	 -S	Append 'singletons' to end of MCL cluster list (excluded by default).
	 -plot Draw a barplot of cluster abundances.

   [dependencies]
	 MCL	(http://www.micans.org/mcl/).
	 R

" }
