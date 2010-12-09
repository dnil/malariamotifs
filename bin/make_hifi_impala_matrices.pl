#!/usr/bin/perl -w

use Bio::Tools::BPlite;
use Bio::Tools::Run::StandAloneBlast;

my $DEBUG = 0;
my $WARNING = 1;
my $LOG = 0;

my $impfile = ""; 
my $sectfile = "";
my $dbname = "";

my $blastpgp_e = 10000;
my $blastpgp_h = 10000;

my $blastpgp_f = 4;
my $blastpgp_F = 'F';
my $blastpgp_b = 1000;
my $blastpgp_v = 1000;

# PLEASE REMEMBER to use dthresh if you want MOTIF SKEW HYPERGEO PVAL THRESHOLD
my $use_dthresh = 0;
my $use_athresh = 0;

while (my $arg = shift @ARGV) {
    if ($arg eq "--db") {
	$arg = shift @ARGV;
	$dbname = $arg;
	$DEBUG && print "DB basename $dbname\n";
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
    elsif ($arg eq "--log") {
	$LOG = 1;
	$DEBUG && print "Using threshold for d-values\n";
    }
    elsif ( $arg eq "--debug" ) {
	$DEBUG = 1;
    }
}

if ($sectfile eq "" or $dbname eq "") {
    print "Missing filename!\n";
    exit;
}

my $hyp_cutoff_threshold = 0.2;

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

$LOG && print "blastpgp_e=$blastpgp_e\n";
$LOG && print "blastpgp_h=$blastpgp_h\n";
$LOG && print "blastpgp_f=$blastpgp_f\n";
$LOG && print "blastpgp_F=$blastpgp_F\n";
$LOG && print "blastpgp_b=$blastpgp_b\n";
$LOG && print "blastpgp_v=$blastpgp_v\n";

# read in motifs as-is from the section file, write as separate fasta files, make blastdbs 
# run blastpgp with one (preferably the original motif seed) as query, checkpoint and make matrices

# clear the sequence entry lists - we're making new ones while rerunning blastpgp..
open SNLIST, ">$dbname.sn"; 
close SNLIST;

open PNLIST, ">$dbname.pn"; 
close PNLIST;

my $ceq_file = "";
my $seed_found = 0;
my $seed_name = "";
my $seed_seq = "";

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

		# this needs the new hypergeo-values to be terribly useful
#		if( $use_dthresh == 1 && $hyperp{$current_motif} >= $hyp_cutoff_threshold ) {		    
#
#		} else {
#		    push @useful_motifs, $current_motif;
#		}

		# close seed db-file
		close FASTA;

		# format database
		system("formatdb -i $local_db_name -p T");

		if($seed_found == 1) {
		    # good, we have a seed!
		} else {
		    print "ERROR: the seed for motif $current_motif was not found among the entries.\n";
		    print "This is potentially ok, but you'll probably want to think about it a while before picking another contig as an interim seed..\n";
		    exit;
		}
		
		run_blastpgp($seed_name, $seed_seq);

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

		    # add ceq-file name to list (kind of slow to use append; the alternative is to keep an open filehandle from caller)
		    open SNLIST, ">>$dbname.sn";
		    print SNLIST $ceq_file."\n";
		    close SNLIST;

		}
	    }
	}
    }
}

# so, that should produce the proper blpgp, chpgp 


# then run the ncbi-toolkit programs to construct the new matrix set..
$readfile_basename = $dbname; # =~ /(.+)\.par\.clustalnames\.nonsingletons/;

system("makemat -P $dbname -d $readfile_basename");
system("copymat -P $dbname -r F");

# END MAIN #

sub run_blastpgp {
    # recieve or construct a Sequence object from the seed to use with the blastreport-factory
    my $seed_name = shift;
    my $seed_seq = shift; 
    
    my $seq = Bio::Seq->new( '-id' => $seed_name, '-seq' => $seed_seq );

    my $local_db_file = $seq->id.".db.fasta"; # ehhrr and else!!
    my $outfile_basename = $local_db_file."_searched_by_".$seed_name;

    # let the motif-file parsing part write the ceq file for us, and update the sn file
    
    # pn is global
    my $chk_file = $outfile_basename.".blpgp.chk"; 

    # add the to be produced chk-file name to list (kind of slow to use append; the alternative is to keep an open filehandle from caller and pass it, or use a global glob)
    # GLOBAL $dbname
    open PNLIST, ">>$dbname.pn";
    print PNLIST $chk_file, "\n";
    close PNLIST;
    
    my @params = ('database' => $local_db_file ,'outfile' => $outfile_basename.".blpgp",
	       'e' => $blastpgp_e, 'h'=>$blastpgp_h, 'f'=>$blastpgp_f, 'F'=>$blastpgp_F, 'b'=>$blastpgp_b, 'v'=>$blastpgp_v,
	       # testade e=1000 med dåliga resultat..
#	       'm'=>6,
	       'C'=>$chk_file,
	       '_READMETHOD'=> 'BPlite'
	       );
# m=6 is a flag to get aln working according to BPlite perldoc

    $DEBUG && print "Blastpgp of ".($seq->id)." (".($seq->seq).") vs $local_db_file..\n";

    my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
    $factory->j(10);
    my $report = $factory->blastpgp($seq);
    
    if( $report->{'REPORT_DONE'} == 1 ) {
	$WARNING && print "OOPS: Report reported done before calling number_of_iterations.\n";
	$WARNING && print "This probably means there were no hits below e-val treshold. Dropping current query.\n";
	return;

	# note that the seed has already been marked status 1..
    }

    my $total_iterations = $report->number_of_iterations;
    my $last_iteration = $report->round($total_iterations);

#    disble the checking for a development test
    
}    


#     my @markupqueue;
#     while (my $sbjct = $last_iteration->nextSbjct) {
	
# 	$DEBUG && print "Iterated blast hit ", $sbjct->name,".\n";
	
# 	while (my $hsp = $sbjct->nextHSP) {
# 	    if($hsp->positive() >= $hsp_pos_cutoff) {
# 		$DEBUG && printf "Motif member %10s  %3d - %3d : %s\n", $hsp->hit->seq_id, $hsp->hit->start, $hsp->hit->end, $hsp->sbjctSeq;
# 		my $hsp_hit_start = $hsp->hit->start;
# 		my $hsp_hit_end = $hsp->hit->end;

# 		# checking the motifs...

# 		push @markupqueue, [$hsp->hit->seq_id, $hsp_hit_start, $hsp_hit_end, 2];

#                 #
# 		# DEVOLP_DEBUG
# 		my $debug_motif_seq_id = $hsp->hit->seq_id;
# 		$debug_motif_seq_id =~ s/\s+$//;
	       
# 		$DEVELOP_DEBUG && print "Retrieving sequence $debug_motif_seq_id.\n";
# 		my ($debug_motif_seq) = ($aln->each_seq_with_id($debug_motif_seq_id));
# 		my $debug_motif = $debug_motif_seq->seq();
# 		$debug_motif =~ s/-//g;
# 		my $debug_motif_substr = substr( $debug_motif, $hsp_hit_start-1, $hsp_hit_end - $hsp_hit_start + 1 );

# 		$DEVELOP_DEBUG && print "Strict section: offset results in motif str ", $debug_motif_substr, ".\n";
# 	    } else {
# 		$DEBUG && print "Dropping HSP since positives ",$hsp->positive," less than $hsp_pos_cutoff (",$hsp->sbjctSeq,") .\n";
# 		# note that this is very powerless, since later motif identification will anyway use all hsps in for the matrix from the .chk-file!?! Better be open about it?
# 		# still useful if all but one have too many mismatches.. 
# 	    }
# 	}
#     }
#     if (@markupqueue > 0) { # the seed is already saved, so the first hit is the second member
# 	# mark subject hit!
# 	while (my $markme = shift @markupqueue) {
# 	    my $save = markup_domain(@$markme);
# 	    if($save == 1) {
# 		# RECHECK - is this correct?
# 		save_motif($seq_name, @$markme);
# 		$DEBUG && print "Saving motif $seq_name ",join(" ",@$markme)," from run_blastpgp for seed ".$$markme[0]."\n";
# 	    }
# 	}
#     } else {
# 	# drop this alignment/motif
# 	$WARNING && print "Number of subjects hit < 1, so dropping current query.\n";
# 	# note that the seed has already been marked status 1.. so don't print it..	
#     }
#     if($strict_section == 1) {
# 	system("rm $local_db_file.phr $local_db_file.pin $local_db_file.psq") unless $DEBUG;
#     }




