#!/usr/bin/perl
use warnings;
use English;
use Tie::File;

#adapted from MLGO (http://www.geneorder.org/).
#this should be run by executing the wrapper script: 'mauve.matrix2adj.pl'

($ARGV[0] and $ARGV[1]) or &HELP_MESSAGE;

my @alphabet;
#my $INPUT_FILE;
my $INPUT_FILE = $ARGV[0];
#my $OUT_FILE;
my $OUT_DIR = $ARGV[1];
my $Dataset;
my $Adjs;
my @Adj;
my @datasetNames = (); #to keep the adjacencies in order

system("mkdir $OUT_DIR");

@alphabet = qw / 0 1 2 3 4 5 6 7 8 9 A B C D E F G H I J K L M N O P Q R S T U V /;

$Dataset = ();
$Adjs = ();
@datasetNames = ();

readFile();
findAllAdjs();

sub findAllAdjs {
    my @concatenations = ();
    my $ind = 0;
    my %hash=();
    foreach my $geName (@datasetNames) {
        my %binaryList = ();
        my@mwadj=();	#for output table of adjacencies.
        foreach my $chroms (@{ $Dataset->{$geName}->[0] }) {
            my @extent = map { $_ > 0 ? (2 * abs($_) - 2, 2 * abs($_)- 1) : (2 * abs($_) - 1, 2 * abs($_) - 2) } @{$chroms};	#transforms all blocks into ordered 2-number 'coordinates'


           #MRW#
            my$mcount=0;
            foreach my$mw (@{$chroms}){
            	my$coords=$extent[$mcount]."-".$extent[$mcount+1];
            	unless(exists( $hash{abs($mw)} )) {
            		$hash{abs($mw)} = $coords;	#matches block number to ordered 2-number 'coordinates'
            	}
            	$mcount+=2;
            }
            #####

            unshift  @extent, -1;
            push     @extent, -1;

            foreach ((1 .. @extent / 2)){
                my $a_b = "$extent[2 * $_ - 1]".",$extent[2 * $_ - 2]"; #these are the adjacenies!
                my $b_a = "$extent[2 * $_ - 2]".",$extent[2 * $_ - 1]";

                if ( exists $Adjs->{$a_b} ){
                	$binaryList{$Adjs->{$a_b}} = 1;	#binary list of adjaceny presence?
                	next;
                }
                if ( exists $Adjs->{$b_a} ){
                	$binaryList{$Adjs->{$b_a}} = 1;
                	next;
                }
                #else:
                $binaryList{$ind} = 1;	#index of where each adj is listed in @Adj
                $Adjs->{$b_a} = $ind;
                push @Adj, $b_a; #seperate list to keep the adjs in order. ALL observed adjacencies!

                $ind ++ ;
            }

            $Dataset->{$geName}->[1] = \%binaryList;


	#### MRW -- save table of adjacencies ####
    	    my$outname = $geName."_adj.txt";
    	    open ADJOUT, ">./$OUT_DIR/$outname";

    	    foreach my$jout (sort keys$$Dataset{$geName}[1]){
    	    	my@sout = split(",",$Adj[$jout]);
    	    	print ADJOUT "Adjacency_".join(",",sort { $a <=> $b }@sout)."\t1\n";
    	    }
    	    close ADJOUT;
    	#### -------------------------------- ####


        }

    #### MRW --- save table of block 'coordinates' ####

    open COOR, ">./$OUT_DIR/block-coords.txt";
    foreach my$hout (sort { $a <=> $b }(keys%hash)){
    	print COOR $hout ."\t". $hash{$hout}. "\n";
    }
    close COOR;
    #### ----------------------------------------- ####


    }


}


sub readFile {
    my $geName;
    my @genes = ();
    open DATASET, "<", $INPUT_FILE || die "can't open file";
    while (<DATASET>) {
        chomp;
        if (/^>(.+)$/){
            $geName = $1;
            next;
        }
        $_ =~ s/\$//g;
        push @{ $Dataset->{$geName}->[0] }, [split /\s/, $_];	# hashref->{key}->[index] = array ?? An array of hashes??
        push @datasetNames, $geName;
    }
    close DATASET;
}

sub HELP_MESSAGE { die "
.Description:
   ATTENTION $0 is not meant to be executed directly, but should be run using the wrapper 'mauve.matrix2adj.pl'

" }





1;
__END__
