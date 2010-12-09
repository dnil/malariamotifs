#!/usr/bin/perl -w

use Bio::Tools::BPlite;

my $DEBUG = 0;

my $impfile = ""; 
my $sectfile = "";
my $oldmotifdir = "";
my $readfile_basename = "";

# PLEASE REMEMBER to use dthresh if you want MOTIF SKEW HYPERGEO PVAL THRESHOLD
my $use_dthresh = 0;

# my $hsp_cutoff_threshold = 50; # used to be 0.001
my $hsp_cutoff_threshold = 680; # does become less important with a subregion filter, but must nevertheless be set? 
my $hyp_cutoff_threshold = 0.09;

while (my $arg = shift @ARGV) {
    if ($arg eq "--db") {
	$arg = shift @ARGV;
	$dbname = $arg;
	$DEBUG && print "DB basename $dbname\n";
    }
    elsif ($arg eq "--basename") {
	$arg = shift @ARGV;
	$readfile_basename = $arg;
	$DEBUG && print "Readfile basename $readfile_basename\n";
    }    
    elsif ($arg eq "--oldmotifs") {
	$arg = shift @ARGV; 
	$oldmotifdir = $arg;
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

if ($readfile_basename eq "") {
    ($readfile_basename) = $dbname =~ /(.+)\.par\.clustalnames\.nonsingletons/;   
}

if ($sectfile eq "" or $dbname eq "" or $oldmotifdir eq "") {
    print "Missing filename!\n";
    exit;
}

open SECTION, $sectfile || exit;
# open(IMPALA, $impfile) || exit;

my $n_motifs;

my $in_motifs = 0;
my $in_seed = 0; 

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

		if( $use_dthresh == 1 && $hyperp{$current_motif} >= $hyp_cutoff_threshold ) {		    

		} else {
		    push @useful_motifs, $current_motif;
		}

		$in_seed = 0;
	    }
	}
    }
}


open SNLIST, ">$dbname.sn";
open PNLIST, ">$dbname.pn";

while ( $motif = shift @useful_motifs ) {

    my $seq_name = $motif;
    my $outfile_basename = $dbname."_searched_by_".$seq_name;
    
    $ceq_file = $outfile_basename.".ceq";    
    system("cp $oldmotifdir/$ceq_file .");

#    open CEQ, ">$ceq_file";
#    print CEQ ">$seq_name\n";
#    print CEQ $seq->seq(),"\n"; #and the sequence was..?!?
#    close CEQ;

    # add ceq-file name to list 
    print SNLIST $ceq_file."\n";

    my $chk_file = $outfile_basename.".blpgp.chk"; 
    system("cp $oldmotifdir/$chk_file .");

#    system("cp $oldmotifdir/$outfile_basename.blpgp.mtx .");

    # add the to be produced chk-file name to list 
    print PNLIST $chk_file, "\n";

}

close SNLIST;
close PNLIST;

# then run the ncbi-toolkit programs to construct the new matrix set..
#($readfile_basename) = $dbname =~ /(.+)\.par\.clustalnames\.nonsingletons/;

system("makemat -P $dbname -d $readfile_basename");
system("copymat -P $dbname -r F");

#
# Usage example:
#
# daniel@ichiban clean_motifs_test]$ ../../malariamotifs/bin/make_new_impala_matrices_for_scoring_scheme.pl --dthresh --section ../malaria.040526.screen.qp.seq.namefix.loq.par.clustalnames.nonsingletons.pep.5dominants.aln.opi.section --db malaria.040526.screen.qp.seq.namefix.loq.par.clustalnames.nonsingletons.5dominants.aln.db.fasta --oldmotifs ../malaria.040526.screen.qp.seq.namefix.loq.par.clustalnames.nonsingletons.5dominants.motifs                    
# 
