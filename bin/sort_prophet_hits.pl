#!/usr/bin/perl -w
# sort_prophet_hits <prophet_file_name> [lowest score to print] 

my @clargs = @ARGV;
my $DEBUG = 0;
my $WARNING = 1;

my $infile = "";
my $topn = 0;
my $threshold = 0;

while (my $arg = shift @clargs) {

    if ($arg eq "--in") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no filename given to --in flag.\n";
	} else {
	    $infile = $arg;
	}
    } elsif($arg eq "--threshold") { # id%?
	$arg = shift @clargs;
	
	if($arg eq "" or $arg <0)  {
	    $WARNING && print STDERR "WARNING: incorrect argument given to --threshold flag.\n";
	} else {
	    $threshold = $arg;
	}
    } elsif($arg eq "--topn") { # id%?
	$arg = shift @clargs;
	
	if($arg eq "" or $arg <= 0)  {
	    $WARNING && print STDERR "WARNING: incorrect argument given to --topn flag.\n";
	} else {
	    $topn = $arg;
	}
    } elsif($arg eq "--debug" or $arg eq "-D") {
	$DEBUG = 1;
	$WARNING && print STDERR "DEBUG mode enabled.\n";
    } elsif($arg eq "--silent" or $arg eq "-q") {
	$WARNING = 0;
    } else {
	$WARNING && print STDERR "WARNING: Unrecognised command line parameter $arg.\n";
	usage();
    }

}

if( $infile eq "") {
    usage();
}

open PROPHET, "<$infile";

my (@hit_entry, @hit_score);
my $in_entry = 0;
my $score = 0;
while($row = <PROPHET>) {
    
    if ($row=~/^Local: Consensus vs .+/) {
#	$current_hit_contig = $1;
# $DEBUG && print $current_hit_contig,"\n";
	if( $in_entry == 1 ) {
	    # save last entry
	    push @hit_entry, $current_entry;	   
	    push @hit_score, $score;
	} else {
	    $in_entry = 1;
	}
	$current_entry = $row;
    } elsif($in_entry == 1) {
	$current_entry .= $row;
	if($row=~/^Score: ([\d\.]+)/) {
	    $score = $1;
	} 
    }    
}
# ..and save the last one..
if($in_entry == 1) {
    push @hit_entry, $current_entry;	   
    push @hit_score, $score;	   
}

my @index = 0..(@hit_entry-1);

my @sorted_index = sort { $hit_score[$b] <=> $hit_score[$a] } @index;

for my $i (@sorted_index) {
    if($hit_score[$i] >= $threshold) {
	push @sorted_hit_entries, $hit_entry[$i];
    }
}

if ($topn == 0) { 
    print join("\n",@sorted_hit_entries);
} else {
    print join("\n",@sorted_hit_entries[0..($topn-1)]);
}

sub usage {
    print "USAGE: sort_prophet_hits.pl\t<--in prophet_file>\n";
    print "\t\t\t\t[--threshold f]\t0<=f\t0\n";
    print "\t\t\t\t[--topn i]\t0<i\tshow all\n";
    exit(1);
}
