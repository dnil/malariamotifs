#!/usr/bin/perl -w

use Bio::Tools::BPlite;

my $impfile = ""; 
my $sectfile = "";

my $DEBUG = 0;

my $use_dthresh = 0;
my $hyp_cutoff_threshold = 0.09;

while (my $arg = shift @ARGV) {
    if ($arg eq "--imp") {
	$arg = shift @ARGV;
	$impfile = $arg;
	$DEBUG && print "Opening $impfile\n";
    }
    elsif ($arg eq "--section") {
	$arg = shift @ARGV;
	$sectfile = $arg;
	$DEBUG && print "Opening $sectfile\n";
    }
    elsif ($arg eq "--dthresh") {
	$use_dthresh = 1;
	$DEBUG && print "Using threshold for d-values\n";
    }
    elsif ($arg eq "-h") {
	$use_dthresh = 1;
	$DEBUG && print "Using threshold for d-values\n";
	$arg = shift @ARGV;
	$hyp_cutoff_threshold = $arg;
    }

    elsif ( $arg eq "--debug" ) {
	$DEBUG = 1;
    }
}

if ($sectfile eq "" or $impfile eq "") {
    print "Missing filename!\n";
    exit;
}

# my $hsp_cutoff_threshold = 50; # used to be 0.001
my $hsp_cutoff_threshold = 680; # used to be 0.001 

open SECTION, $sectfile || exit;
open(IMPALA, $impfile) || exit;

my $n_motifs;

my $in_motifs = 0;
my $in_seed = 0;    # the in_seed is overkill so far, but might come in handy at one point

my %se_count;
my %co_count;
my %se_fract;
my %hyperp;
my %skewdir;
my %nmembers;

my $current_motif;

while (my $row = <SECTION>) {
    chomp $row;
    if ($row !~ /\w/) {
    } elsif($in_motifs == 0) {
	if( $row =~ m/^Found (\d+) motifs\./) {
	    $n_motifs = $1;

	    $in_motifs = 1;

	    $DEBUG && print $row,"\n";
	}
    } else {	
	if($in_seed == 0) {
	    if ($row =~ m/Seed\s+(\S+)\s+-\s+(\d+)\s+members/ ) {
		$in_seed = 1;

		$current_motif = $1;
		$nmembers{$current_motif} = $2;
		$DEBUG && print $row,"\n";
	    }
	} else {
	    # in Seed 	    
	    if($row =~ m/SE\s+(\d+)\s+CO\s+(\d+)\s+SEfraction\s+([\d\.]+)\s+p\s+([\d\.]+)\s+([+-])/) {
		$se_count{$current_motif} = $1;
		$co_count{$current_motif} = $2;
		$se_fract{$current_motif} = $3;
		$hyperp{$current_motif} = $4;
		$skewdir{$current_motif} = $5;

		$DEBUG && print $hyperp{$current_motif}, "\t", $se_fract{$current_motif}, "(",$skewdir{$current_motif},")\n";

		$in_seed = 0;
	    }
	}
    }
}

my $impala_report = new Bio::Tools::BPlite(-fh=>\*IMPALA);

#    print $report->query; 
my $logsumscore = 0;
my $some_other_score = 0;

my $se_score = 0;
my $mild_score = 0; 

while( my $sbjct = $impala_report->nextSbjct ) {

    while( $hsp = $sbjct->nextHSP ) {

	my ($hspmotif) =($hsp->hit->seq_id);
	$hspmotif =~ s/\s+$//;
	
	if(defined($hyperp{$hspmotif})) {
	    	    print sprintf("%-20s\t",$hspmotif),$hsp->query->seq_id, "\t", $hsp->query->start, "\t", $hsp->query->end , "\t", $hsp->EXP, "\t", 
	    $hyperp{$hspmotif}, "\t", $se_fract{$hspmotif}, "(",$skewdir{$hspmotif},")\t",$nmembers{$hspmotif},"\n";

	    if($hsp->EXP < $hsp_cutoff_threshold) {
		
		if($use_dthresh == 1 && $hyperp{$hspmotif} >= $hyp_cutoff_threshold ) {
		    # hypergeo p-value cutoff
		    next;
		}
		
		$blastpval = 1 - exp(- $hsp->EXP );

		if($skewdir{$hspmotif} eq "+") {
		    $some_other_score += log($blastpval)+log($hyperp{$hspmotif});
		    $se_score += log($blastpval)+log($hyperp{$hspmotif});
		} else {
		    $some_other_score -= log($blastpval)+log($hyperp{$hspmotif});
		    $mild_score -= log($blastpval)+log($hyperp{$hspmotif});
		}
	    }
	} else {
	    $DEBUG && print "No defined motif $hspmotif -- presume dropped... Eliminate those during section build!\n";
	    print $hspmotif,"\t\t", $hsp->query->seq_id, "\t", $hsp->query->start, "\t", $hsp->query->end , "\t", $hsp->EXP, "\n";	    
	}
    }
}

if( $use_dthresh == 1) {
    print "S = $some_other_score, M = $mild_score, SE = $se_score [hsp->EXP cutoff $hsp_cutoff_threshold, motif skew hypergeo cutoff $hyp_cutoff_threshold]\n";
} else {
    print "S = $some_other_score, M = $mild_score, SE = $se_score [hsp->EXP cutoff $hsp_cutoff_threshold]\n";
}
