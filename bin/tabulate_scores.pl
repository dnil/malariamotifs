#!/usr/bin/perl -w

while($row = <STDIN>) {

    if($row =~m/cutoff/) {
	
	$row =~ m/^(.+)\.fasta.+\:S = ([-]*\d+\.*\d*), M = ([-]*\d+\.*\d*), SE = ([-]*\d+\.*\d*) /;
	
	print $1,"\t",$2, "\t", $3, "\t", $4,"\n";
    }
}
