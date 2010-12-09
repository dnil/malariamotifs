#!/usr/bin/perl -w

use Bio::Tools::BPlite;
use Bio::Tools::Run::StandAloneBlast;

my $DEBUG = 0;
my $WARNING = 1;
my $LOG = 0;

my $impfile = ""; 
my $sectfile = "";
my $dbname = "";

# PLEASE REMEMBER to use dthresh if you want MOTIF SKEW HYPERGEO PVAL THRESHOLD
#my $use_dthresh = 0;
#my $use_athresh = 0;

while (my $arg = shift @ARGV) {
#    if ($arg eq "--db") {
#	$arg = shift @ARGV;
#	$dbname = $arg;
#	$DEBUG && print "DB basename $dbname\n";
#   }
    if ($arg eq "--section") {
	$arg = shift @ARGV;
	$sectfile = $arg;
	$DEBUG && print "Opening $sectfile\n";
    }
#    elsif ($arg eq "--dthresh") {
#	$use_dthresh = 1;
#	$DEBUG && print "Using threshold for d-values\n";
#    }
    elsif ($arg eq "--log") {
	$LOG = 1;
    }
    elsif ( $arg eq "--debug" ) {
	$DEBUG = 1;
    }
}

if ($sectfile eq "") {
    print "Missing filename!\n";
    exit;
}

# my $hyp_cutoff_threshold = 0.2;

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

#$LOG && print "blastpgp_e=$blastpgp_e\n";
#$LOG && print "blastpgp_h=$blastpgp_h\n";
#$LOG && print "blastpgp_f=$blastpgp_f\n";
#$LOG && print "blastpgp_F=$blastpgp_F\n";
#$LOG && print "blastpgp_b=$blastpgp_b\n";
#$LOG && print "blastpgp_v=$blastpgp_v\n";

# read in motifs as-is from the section file, write as separate fasta files, make blastdbs 
# run blastpgp with one (preferably the original motif seed) as query, checkpoint and make matrices

# clear the sequence entry lists - we're making new ones while rerunning blastpgp..
#open SNLIST, ">$dbname.sn"; 
#close SNLIST;

#open PNLIST, ">$dbname.pn"; 
#close PNLIST;

my $ceq_file = "";
my $seed_found = 0;
my $seed_name = "";
my $seed_seq = "";

open MOTIF_FASTA_LIST, ">$sectfile.motif_fasta_list";

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
		$seed_name = $current_motif; # retains _start_stop coordinates in name

		$nmembers{$current_motif} = $2;

		if( $current_motif =~ /^([^_]+)_.+/ ) { 
		    $current_motif = $1;
		}

		$DEBUG && print $row,"\n";

		# open a new seed db-file 
		$local_db_name = "$seed_name.db.fasta"; 

		open FASTA, ">$local_db_name";
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

		# close seed db-file
		close FASTA;

		# add this fasta file to a list of dumped sequences
		print MOTIF_FASTA_LIST "$local_db_name\n"; 

		# format database
#		system("formatdb -i $local_db_name -p T");

		if($seed_found == 1) {
		    # good, we have a seed!
		} else {
		    print "ERROR: the seed for motif $current_motif was not found among the entries.\n";
		    print "This is potentially ok, but you'll probably want to think about it a while before picking another contig as an interim seed..\n";
		    exit;
		}
		
#		run_blastpgp($seed_name, $seed_seq);

#		check_that_we_have_all_members_in_new_hitlist();

		$in_seed = 0;
		$seed_found = 0;
	    } elsif ( $row =~ /^\s*(\S+)\s+(\d+)\s+(\d+)\s+(\w+)\s+(\S+)/ ) {
		my $motif_member = $1;
		my $motif_member_start = $2;
		my $motif_member_end = $3;
		my $motif_member_seq = $4;
		my $motif_member_dom = $6;
		# print something nice to a fasta file or to a string for later output?

		# to dbfile-fasta...
		print FASTA ">".$motif_member."\n"; # "\t".$motif_member_start."\t".$motif_member_end."\n";
		print FASTA $motif_member_seq,"\n";

		if ( $current_motif eq $motif_member ) {
		    $DEBUG && print "DEBUG: found $motif_member to be seed $current_motif with sequence $motif_member_seq.\n";

		    # this is the seed!
		    $seed_found = 1;		    

		    $seed_seq = $motif_member_seq;
		    
		    # write seed to a ceq file
		    my $seed_name =  $motif_member."_".$motif_member_start."_".$motif_member_end;
		    $ceq_file = "$local_db_name.ceq"; # "_searched_by_$seed_name.ceq";

		    open CEQ, ">$ceq_file";
		    print CEQ ">".$seed_name."\n";
		    print CEQ "$motif_member_seq\n";
		    close CEQ;

		}
	    }
	}
    }
}

close MOTIF_FASTA_LIST;
