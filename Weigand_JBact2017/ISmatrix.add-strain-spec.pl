#!/usr/bin/perl -w
use strict;
use Getopt::Long;

&GetOptions(	'ogs=s' => \my$ogs,		# input table of ogs
		'spec=s' => \my$spec,		# output list of strain-specific insertions
		'out=s' => \my$ogs_out, # new output table of ogs.
		'q' => \my$q, # run quietly.
		'h' => \my$h,
		'nondir=s' => \my$nondir);		# directory of unmatched insertions

($ogs and $spec and $nondir and $ogs_out) or &HELP_MESSAGE;
if($h){ &HELP_MESSAGE};

#compile count of unmatched IS insertions in each genome.
open SPEC, ">$spec";
my%strainspec=();
foreach my$nomatch (<"$nondir/*txt">){
        #print $i ."\n";
        open IN, "$nomatch";
        while(my$noj = <IN>){
                chomp$noj;
                #print $j . "\n";
                $strainspec{$noj}=1;
        }
}
#print ("\n",(keys%hash))."\n";
foreach my$no (sort(keys%strainspec)){
        my$no_count = qx(grep -c \"$no\" $ogs);
				chomp$no_count;
				if(scalar($no_count) == 0){
					print SPEC $no."\t".$no_count."\n";
				}
}
close SPEC;

#add remaining unmatched insertions to existing matrix
my%hash=();
my%hout=();
open LIST, "$spec";
while(my$i = <LIST>){
	chomp$i;
	my@si=split("\t",$i);
	my$pdl = substr($i,0,4);
	my($c1)= $si[0] =~ /(?<=\d\.)(\d+)(?=-\d)/;
	my($c2)= $si[0] =~ /(?<=\d-)(\d+)(?=$)/;
	#print join("--",($pdl,$c1,$c2,$si[0]))."\n";
	$hash{$pdl}{$c1}=$si[0];
	$hash{$pdl}{$c2}=$si[0];
	$hout{$si[0]}=1;
}
close LIST;

my@cc=();
open MATRIX, "$ogs";
open OUT, ">$ogs_out";
unless($q){print "Found matches:\n"};
my@header=();
while(my$j = <MATRIX>){
	unless($j =~ /IS481|IS1663|IS1002/){
		print OUT $j; #header
		chomp$j;
		@header = split("\t",$j);
	}else{
		chomp$j;
		my@sj=split("\t",$j);
		my$count=0;
		foreach my$k (@sj){
			if((length($k) < 2) || ($k !~ /IS481/)){
				#skip "-"

			}else{
				my$id=substr($k,0,4);
				if(exists($hash{$id})){
					my@cc = keys$hash{$id};

					my@skk = split(",",$k);
					my@ccc=();
					foreach my$skkk (@skk){
						my($cc1)= $k =~ /(?<=\d\.)(\d+)(?=-\d)/;
						my($cc2)= $k =~ /(?<=\d-)(\d+)(?=$)/;
						push(@ccc,$cc1);
						push(@ccc,$cc2);
					}
					foreach my$c5 (@cc){ #strain-spec list
						foreach my$c6 (@ccc){ #matched in matrix
							if((abs($c5 - $c6) <= 10)){

								my($cout1) = $hash{$id}{$c5} =~ /(?<=\d\.)(\d+)(?=-\d)/;
								my($cout2) = $hash{$id}{$c5} =~ /(?<=\d-)(\d+)(?=$)/;

								if(($c5 <= $c6) && ($cout1 > $cout2)){
									unless($q){print "\t". $k.",".$hash{$id}{$c5}."\n"};
									$sj[$count] = $k.",".$hash{$id}{$c5};
								}elsif(($c5 <= $c6) && ($cout1 < $cout2)){
									unless($q){print "\t". $hash{$id}{$c5}.",".$k."\n"};
									$sj[$count]= $hash{$id}{$c5}.",".$k;
								}elsif(($c5 >= $c6) && ($cout1 > $cout2)){
									unless($q){print "\t". $hash{$id}{$c5}.",".$k."\n"};
									$sj[$count]=$hash{$id}{$c5}.",".$k;
								}else{
									#print $k."\t|\t".$hash{$id}{$c5}."\n";
								}
								#remove from list any that match as 'duplicates'
								delete$hout{ $hash{$id}{$c5} };
								last;
							}
						}
					}
				}
			}
			$count++;
		}
		print OUT join("\t",@sj)."\n"
	}
}

unless($q){print "Unmatched, strain-specific insertions:\n"};
my@remain = sort(keys%hout);
for (my$h=0; $h < @header; $h++){
	my@end = ("-") x scalar(@header);
	foreach my$oo (@remain){
		if($oo =~ /$header[$h]/){
			unless($q){print "\t". $header[$h]."--".$oo."\n"};
			$end[$h] = $oo;
			print OUT join("\t",@end)."\n";
		}
	}

}

##########################
sub HELP_MESSAGE { die "
.Description:
   Takes a directory of unmatched IS elements from 'ISmatrix2rbm.pl' and adds potentiall strain-specific IS elements to an existing matrix from 'ogs.rb'.

.Usage: $0 -nondir -spec -ogs -out

[mandatory]
	-nondir	<path>		Directory of unmatched IS insertions.
				Expects containing files to be named like: 'A001_CP000001-A002_CP000002-unmatched.txt'
	-spec	<spec.txt>	List of identified strain specific insertions.
	-ogs	<in.ogs>	Matrix of IS matches from 'ogs.rb'.
	-out	<out.ogs>	New matrix of IS matches with added strain-specific insertions.

[optional]
	-q	Run quietly.
	-h	This helpful message.

WARNING: Written to specifically recognize name formats containing PDL strain-ids and GenBank accession numbers, eg. 'A001_CP000001'.
 
" }
