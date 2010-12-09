#!/usr/bin/perl -w

my $DEBUG = 0;
my $WARNING = 1;

my $fastafile = "";
my $set1_fract = 2/3; 

while (my $arg = shift @ARGV) {

    if ($arg =~ /^-/) {	
	if($arg eq '-D' or $arg eq '--debug') { 
	    $DEBUG = 1;
	} elsif ($arg eq '-p') {
	    $arg = shift @ARGV;
	    $set1_fract = $arg;
	}
    } else {
	$fastafile = $arg;
    }
}

if($fastafile eq "") {
    print "No fasta file to split was given. Nothing to do!\n";
    exit 1;
}

open FASTAFILE, $fastafile or die "Could not open $fastafile.\n";

my $set2_fract = 1 - $set1_fract;

my $printentry = 0;

my %seq; 
my %header;

my $current_nr = -1;

while( my $row = <FASTAFILE> ) {    
    chomp $row;
    if ( $row =~ m/^>/ ) {
	$DEBUG && print STDERR "Checking $row.\n";
	if(($current_nr) = ($row =~m/^>\D+(\d+)/)) {
	   $DEBUG && print $current_nr,"\n";
	    
	    $header{$current_nr} = $row;
	    $seq{$current_nr} = "";

	} else {
	    $WARNING && print "WARNING: Fasta header row did not contain the magic contig number - ignoring $row.\n";
	}
    } else {
	if ($current_nr != -1) {
	    $seq{$current_nr} .= $row;
	    $DEBUG && print "added ", length($row), "bases to sequence $current_nr.\n";
	}
    }
}

my @numbers = keys %header;
my $entries = @numbers;

for (my $i = 0; $i < $entries; $i++) {   
    $pick_next = int(rand($entries - $i));

    push @new_order, $numbers[$pick_next];
    splice @numbers, $pick_next, 1;
}

if (@numbers > 0 ) {
    print "Oops, ",scalar(@numbers)," numbers left in randomisation queue!";
}

my $setfh; 
open PFH, ">$fastafile.p.set";
open QFH, ">$fastafile.q.set";

my $outfh = "PFH";
my $pcount = 0;
my $qcount = 0;

while ( my $next = shift @new_order ) {

    print $outfh $header{$next}, "\n";

    my $linecount=0;
    my @nextseq = split (/ */,$seq{$next});

    $DEBUG && print "DEBUG: length of sequence $next is ",length($seq{$next}), " or ", scalar(@nextseq),"\n";       
    
    while( $nextbase = shift @nextseq ) {
	$linecount++;
	print $outfh $nextbase;
	if($linecount == 60) {
	    print $outfh "\n";
	    $linecount = 0;
	}

    }
    if ($linecount != 0) {
	print $outfh "\n";
    }

    if($pcount/$entries > $set1_fract ) {
	$outfh = "QFH";
	$qcount++;
    } else {
	$pcount++;
    }
}

close PFH;
close QFH;

print STDERR "Order randomised. Wrote $pcount entries to $fastafile.p.set and $qcount to $fastafile.q.set.\n";
