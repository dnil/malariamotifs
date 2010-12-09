#!/usr/bin/perl -w

my %ids; 
while ($row = <STDIN>) {
    chomp $row; 
    ($score) = ($row =~ /S = ([+-\d\.]+)/);
    (@isolates) = ($row =~ /^([^:]+)/);
    while ($i = shift @isolates) { 
	$ids{$i} = $score;
    }
} 

@keys = sort keys %ids; 
while ($k = shift @keys) {
    print $k, "\t", $ids{$k},"\n";
}
