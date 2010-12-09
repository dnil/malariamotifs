#!/usr/bin/perl -w

#use Bio;
use Bio::Tools::BPlite;

my $impfile = ""; 

while (my $arg = shift @ARGV) {
    if ($arg eq "--imp") {
	$arg = shift @ARGV;
	$impfile = $arg;
	print "Opening $impfile\n";
    }
}

open(IMPALA, $impfile);
    
my $impala_report = new Bio::Tools::BPlite(-fh=>\*IMPALA);

#    print $report->query;
while( my $sbjct = $impala_report->nextSbjct ) {
    while( $hsp = $sbjct->nextHSP ) {
	print  $hsp->hit->seq_id, "\t", $hsp->query->seq_id, "\t", $hsp->query->start, "\t", $hsp->query->end , "\t", ,"\n";
    }
}
