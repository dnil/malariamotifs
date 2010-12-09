#!/usr/bin/perl -w

use Bio::Tools::BPlite;

my $impfile = ""; 
my $sectfile = "";

my $DEBUG = 0;

my $use_dthresh = 0;

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
my $hyp_cutoff_threshold = 0.2;

my $tolerated_overlap_with_previous_hits = 3;

open SECTION, $sectfile || exit;
open(IMPALA, $impfile) || exit;

my $n_motifs;

my $in_motifs = 0;
my $in_seed = 0; 

my %se_count;
my %co_count;
my %se_fract;
my %hyperp;
my %skewdir;
my %nmembers;

# sort hits in length order, longest first

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

my $query_len = $impala_report->qlength;
my @hitstatus = split (/ */,'0'x$query_len);

my @hitseed = split (/ */,'.'x$query_len);

while( my $sbjct = $impala_report->nextSbjct ) {

    # initialise hit markup string
    while( $hsp = $sbjct->nextHSP ) {

	my ($hspmotif) =($hsp->hit->seq_id);
	$hspmotif =~ s/\s+$//;
	
	if(defined($hyperp{$hspmotif})) {
	    print "CANDIDATE: ",sprintf("%-20s\t",$hspmotif),$hsp->query->seq_id, "\t", $hsp->query->start, "\t", $hsp->query->end , "\t", $hsp->EXP, "\t",
	                        $hyperp{$hspmotif}, "\t", $se_fract{$hspmotif}, "(",$skewdir{$hspmotif},")\t",$nmembers{$hspmotif},"\n";

	    if($hsp->EXP < $hsp_cutoff_threshold) {
		
		if($use_dthresh == 1 && $hyperp{$hspmotif} >= $hyp_cutoff_threshold ) {
		    # hypergeo p-value cutoff
		    print "DROP ",$hyperp{$hspmotif}," >= $hyp_cutoff_threshold.\n";
		    next;
		}

		# the motif has passed the dthresh. time to check if the hit region has already been visited by 
		# a skewed enough motif.

		$blastpval = 1 - exp(- $hsp->EXP );
		$logsum = log($blastpval)+log($hyperp{$hspmotif});

		if (check_domain_score($hsp->query->start, $hsp->query->end, $logsum) != 1) {
		    # attempt to lock
		    markup_domain($hsp->query->start, $hsp->query->end, $logsum, $hspmotif);
		    
		} else {
		    print "DROP: better hit between ",$hsp->query->start," ",$hsp->query->end," than $logsum already.\n";
		}
	    }
	} else {
	    $DEBUG && print "No defined motif $hspmotif -- presume dropped... Eliminate those during section build!\n";
	    print $hspmotif,"\t\t", $hsp->query->seq_id, "\t", $hsp->query->start, "\t", $hsp->query->end , "\t", $hsp->EXP, "\n";
	}
    }
}

my $current_seed = "";
for ( my $i = 0 ; $i < $query_len; $i++ ) {
    my $hspmotif = $hitseed[$i];

    if ($hspmotif eq ".") {
	$DEBUG && print "Noting at pos ",($i+1),".\n";
	next;
    }

    if( $hspmotif ne $current_seed ) {
	my $logsum = $hitstatus[$i];		       

	$DEBUG && print "Final hit for $hspmotif used with logscore $logsum (skewdir ",$skewdir{$hspmotif},") starting at pos ",($i+1),"\n";

	if($skewdir{$hspmotif} eq "+") {
	    $some_other_score += $logsum;
	    $se_score += $logsum;
	} else {
	    $some_other_score -= $logsum;
	    $mild_score -= $logsum;
	}	

	$current_seed = $hspmotif;
    } else {
	$DEBUG && print "Nothing new to see at pos ",($i+1),"(still ",$hspmotif,")\n";
    }    
}

if( $use_dthresh == 1) {
    print "S = $some_other_score, M = $mild_score, SE = $se_score [hsp->EXP cutoff $hsp_cutoff_threshold, motif skew hypergeo cutoff $hyp_cutoff_threshold]\n";
} else {
    print "S = $some_other_score, M = $mild_score, SE = $se_score [hsp->EXP cutoff $hsp_cutoff_threshold]\n";
}

sub check_domain_score {
    # basically an array, where each element holds a P val
    # if hitstat is set for some positions, check there for best or worst P val to compare current hit with

    my $start = (shift) - 1; #convert to internal 0-based desc.
    my $end = (shift) - 1;
    my $pval = shift;

    my $len = $end - $start + 1;

#    my @hs = @hitstatus[$start, $end];
    
    my $hitstatusoverride_count = 0;

    for (my $i = $start; $i<= $end ; $i++) {
	if ($hitstatus[$i] > $pval) {
	    $hitstatusoverride_count++;
	}
    }

#    if (     $WARNING && print "WARNING: only a part of previous hit was better than current: ".$hsp[$i]." at pos ".($i+1)." compared to current $pval\n";

    # many ways to do the greedy overlay; one is this, where we need a certain least number of matches to overwrite for the entire length of the motif
    # an other, more conservative?, approach would be to not update if there were a certain number of old hit columns within the motif with better pvals

    if ( $hitstatusoverride_count == 0 ) {
	# all old were better than new. case closed.
	$DEBUG && print "Nothing to see: hitstatusoverride_count == 0. No news.\n";
	return(1);
    } elsif ( $hitstatusoverride_count < $tolerated_overlap_with_previous_hits ) {
	$DEBUG && print "Less than $tolerated_overlap_with_previous_hits changes, $hitstatusoverride_count, out of $len to $start - $end hitstatus would occur. Ignoring the new, partially better, hit.\n";	
	return(1);
#    } elsif ( $len - $hitstatusoverride_count ) {
    } else {
	$DEBUG && print "DEBUG Checking status for $start - $end: New is better (there were $hitstatusoverride_count positions of new motif len $len with better pvals than previous hit).\n";
	return(0);
    }
}

sub markup_domain {
    # basically a string of chars
    # 0: untested, unhit
    # 1: tested
    # 2: hit by other element

    my $start = (shift) - 1; #convert to internal 0-based desc.
    my $end = (shift) - 1;
    my $pval = shift;
    my $seedname = shift;

    my $hit_len = $end - $start + 1;

    # test number of overrides that would ensue from 
#    my $hitstatusoverride_count = 0;

#    for(my $i = $start; $i <= $stop; $i++) {
#	if( $hitstatus[$i] > $pval ) {
#	    $hitstatusoverride_count++;	    
#	}	
#    }


    # actual replace

    for(my $i = $start; $i <= $end; $i++) {
	$hitstatus[$i] = $pval;
	$hitseed[$i] = $seedname;
    }
    
    $DEBUG && print "setting status for $start - $end to $pval.\n";

    return 0;
}
