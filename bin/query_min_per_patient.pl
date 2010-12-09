#!/usr/bin/perl -w

my %ids; 
while ($row = <STDIN>) {
    chomp $row; 
    ($score) = ($row =~ /S = ([-+\d\.]+)/);
    (@isolates) = ($row =~ /_(\w{2,3}\d{1,2})\.\d_/g);
    while ($i = shift @isolates) { 
	if(!defined($ids{$i}) || $score < $ids{$i}) {
	    $ids{$i} = $score;
	}
    }
} 

@keys = sort keys %ids; 
while ($k = shift @keys) {
    print $k, "\t", $ids{$k},"\n";
}
