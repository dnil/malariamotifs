#!/usr/bin/perl -w

my $last = ""; 

while(<>) { 
    chomp; 
    @r = split /\s+/; 
    if ($last eq $r[2]) { 
	print "\t\t$r[0] \t($r[1])\t$r[3]\t$r[6]\t$r[5]\t$r[7]\n"; 
    } else { 
	print "$r[2]\t$r[4]\t$r[0] \t($r[1])\t$r[3]\t$r[6]\t$r[5]\t$r[7]\n"; 
	$last = $r[2]; 
    }
}
# > malaria.040526.screen.qp.fasta.ace.fractioncount.sorted_out
