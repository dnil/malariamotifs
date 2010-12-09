#!/usr/bin/perl -w

# v2: doesn't have to count dominance, reads dominance from file
#     so, no dominance exclusion, or blast stuff needed then
#     should really not need assembly or other info about the individual peptides
#     other than a plain alignment and the dominance file

use Bio::AlignIO;
use Bio::SearchIO;
use Bio::SeqIO;

use Bio::Tools::Run::StandAloneBlast;

sub read_dominance_table;

sub find_peak;
sub check_remaining_hits_in_region;
sub run_blastpgp;
sub markup_domain;
sub check_domain_status;
sub get_unchecked_parts;

sub output;
sub usage;

my @clargs = @ARGV;
my $DEBUG = 0;
my $WARNING = 1;

my $DEVELOP_DEBUG = $DEBUG;

my $LOG = 0;

my $SE_UKS_COMPAT = 1;

my $outfile = "";
my $infile = "";

my $dominance_table_name = "";

# apparent blastpgp limitation (Or is that simply due to the "no hits found" bioperl issue?)
my $min_blastp_len = 6;

# default values
#my $min_init_w = 4;
#my $maxw = 40;
#my $ignore_section_treshold = 0;

#my $ultraconserved_threshold = 0.9;
my $ultraconserved_threshold = 0.82;
#my $conserved_threshold = 0.45;
my $conserved_threshold = 0.4;

my $tolerated_motif_overlap = 3; # a tad high for 6 aa motifs? 
# possibly good for additional hits; very annoying for seed picking!

my $hsp_neg_cutoff = 2; # isn't that a bit harsh?

my $count_one_dominant_per_isolate = 0;

my $strict_section = 0;

#v3..
my $status_aware = 1;

# blastpgp run values

my $blastpgp_e = 100;
my $blastpgp_h = 10;

my $blastpgp_f = 4;
my $blastpgp_F = 'F';
my $blastpgp_b = 1000;
my $blastpgp_v = 1000;

while (my $arg = shift @clargs) {

    if ($arg eq "--in") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no filename given to --in flag.\n";
	} else {
	    $infile = $arg;
	}
    } elsif($arg eq "--ultraconserved") { # id%?
	$arg = shift @clargs;
	
	if($arg eq "" or $arg <0 or $arg >1)  {
	    $WARNING && print STDERR "WARNING: incorrect argument given to --ultraconserved flag.\n";
	} else {
	    $ultraconserved_threshold = $arg;
	}
    } elsif ($arg eq "--countmodel") {
	$arg = shift @clargs;
	if ($arg eq "oneperisolate") {
	    $count_one_dominant_per_isolate = 1;
	} else {
	    print STDERR "WARNING: unknown argument given to --countmodel.\n";	    
	}
    } elsif($arg eq "--conserved") { # id%?
	$arg = shift @clargs;
	
	if($arg eq "" or $arg <0 or $arg >1)  {
	    $WARNING && print STDERR "WARNING: incorrect argument given to --conserved flag.\n";
	} else {
	    $conserved_threshold = $arg;
	}
    } elsif($arg eq "--tolerated_overlap") {
	$arg = shift @clargs;
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no argument given to --tolerated_overlap flag.\n";
 	} else {
	    $tolerated_motif_overlap = $arg;
	}
    } elsif($arg eq "--dominancetab" or $arg eq "-t") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no filename given after --dominancetab flag.\n";
	} else {
	    $dominance_table_name = $arg;
	}
    } elsif($arg eq "--blastpgp_negatives") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no rows to delete given to --blastpgp_negatives flag.\n";
	} else {
	    $hsp_neg_cutoff = $arg;
	}
    } elsif($arg eq "--strictsection" or $arg eq "-s") {
	$strict_section = 1;	
	$WARNING && print STDERR "Strict section mode enabled.\n";
    } elsif($arg eq "--debug" or $arg eq "-D") {
	$DEBUG = 1;
	$WARNING && print STDERR "DEBUG mode enabled.\n";
    } elsif($arg eq "--log" or $arg eq "-l") {
	$LOG = 1;
	$WARNING && print STDERR "Log mode enabled.\n";
    } elsif($arg eq "--silent" or $arg eq "-q") {
	$WARNING = 0;
    } else {
	$WARNING && print STDERR "WARNING: Unrecognised command line parameter $arg.\n";
	usage();
    }
}

if( $infile eq "" ) {
    $WARNING && print STDERR "FATAL: Missing required parameter infile.\n";
    usage();
}

my @section_log;

# LOG execution parameters
output "$0 ".(join(" ", @ARGV))."\n";
output "Run at ".(scalar(localtime))." on ".($ENV{HOSTNAME})." in ".($ENV{PWD})." by ".($ENV{USER})."\n";

output "min_blastp_len=$min_blastp_len\n";

output "countmodel oneperisolate=$count_one_dominant_per_isolate\n";
output "strict_section=$strict_section\n";
output "ultraconserved_threshold=$ultraconserved_threshold\n";
output "conserved_threshold=$conserved_threshold\n";

output "tolerated_motif_overlap=$tolerated_motif_overlap\n";
output "hsp_neg_cutoff=$hsp_neg_cutoff\n";

output "blastpgp_e=$blastpgp_e\n";
output "blastpgp_h=$blastpgp_h\n";
output "blastpgp_f=$blastpgp_f\n";
output "blastpgp_F=$blastpgp_F\n";
output "blastpgp_b=$blastpgp_b\n";
output "blastpgp_v=$blastpgp_v\n";

my @exclude_list;
my $exclude = 0;

$ain = Bio::AlignIO->new('-file'=>$infile, '-format'=>'clustalw');
# $aout = Bio::AlignIO->new('-file'=>">$outfile", '-format'=>'clustalw');

my $aln = $ain->next_aln();

my $aln_start = 1; 
my $aln_end = $aln->length;

my $aln_height = $aln->no_sequences;

# read in the dominance table file, if given 

if ($dominance_table_name eq "")  {
    print STDERR "FATAL: no dominance table given.\n";
    exit;
}

my %fractcount;
my %isolatetype_dominants;

open DT, "<$dominance_table_name";
read_dominance_table(\*DT); 
close DT;

print STDERR "\nFrom fraction-count:";
foreach my $isolatetype ( keys %isolatetype_dominants ) {
    print STDERR " ",uc($isolatetype)," dominants: ", $isolatetype_dominants{$isolatetype};
}
print STDERR ".\n\n";

#output "\nFrom fraction-count: CO dominants: $co_dominants SE dominants: $se_dominants.\n\n";

my $db_file = $infile.".db.fasta";
    
# prepare for later psiblast by making a blast-db (unless strict, which makes separate local db)

if( $strict_section == 0 ) {
    # normal, non strict
    $db_temp_seqio = Bio::SeqIO->new('-file'=>">$db_file", '-format'=>'Fasta');
    foreach my $seq ($aln->each_seq) {
	my $seq_str = $seq->seq();
	$seq_str =~ s/-//g;
	$outseq = Bio::Seq->new( '-id' => $seq->id, '-seq' => $seq_str );
	$db_temp_seqio->write_seq($outseq);
    }

    system("formatdb -i $db_file -p T");
}

# clear the seed and matrix out name list files
open CLEAR, ">".$db_file.".sn";
close CLEAR;

open CLEAR, ">".$db_file.".pn";
close CLEAR;

# debug log

my @maxcount;
my @gap;
my @peak;
my @idp;
my @lidp;
my @res;

for (my $i = $aln_start; $i <= $aln_end; $i++) {
    
    my %dominance_by_res;

    my $total_count = 0;
    my %count;

    foreach my $seq ($aln->each_seq) {
	$res = substr($seq->seq,$i-1,1); # aln coords are 1-based inclusive
	$count{$res}++;
	$total_count++;	
	if(defined($fractcount{$seq->id}) && $fractcount{$seq->id} ne "") {
	    $dominance_by_res{$res} .= " ".$fractcount{$seq->id};
	}
    }

    my $max_count = 0;
    my $max_res = '?';
    my $gap_count = 0;

    foreach $res (keys %count) {
	if($res eq '-') {
	    $gap_count = $count{$res};
	} elsif($count{$res} > $max_count) {
	    $max_count = $count{$res};
	    $max_res = $res;
	}
	if( defined($dominance_by_res{$res}) ) {
	    output sprintf("Column: %3d  Res: %s  Count: %3d  Dom: %s\n", $i, $res, $count{$res}, $dominance_by_res{$res});
	} else {
	    output sprintf("Column: %3d  Res: %s  Count: %3d\n", $i, $res, $count{$res});
	}
    }

    push @res, $max_res;

#    if($gap_count > $max_count) {
    if( $gap_count > ($total_count - $gap_count) ) { # gap if majority is gap.
	push @gap, 1;
    } else {
	push @gap, 0;
    }

    $id_percent = $max_count / ($total_count - $gap_count);
#    $id_height = $max_count;
# oops on all gap.. please fix!

    if( $max_count / $aln_height > $ultraconserved_threshold ) {
	push @peak, 1;
    } else {
	push @peak, 0;
    }

    if($id_percent > $conserved_threshold) {
	push @idp, 1;
    } else {
	push @idp,0;
    }
    
    # test a third level - probably not all that useful
    if($id_percent > .15) {
	push @lidp, 1;
    } else {
	push @lidp,0;
    }

    output sprintf("\tGap: %1d  Gaps: %3d  Max count: %3d\n",$gap[$i-1],$gap_count, $max_count);
    output sprintf("\tid: %.3f\n",$id_percent,$max_count);
    output "\n";
}

# smooth peak, idp and gap

my $peak = join('',@peak);
my $gap =  join('',@gap);
my $idp = join('',@idp);
my $lidp = join('',@lidp);

$peak =~ s/10(?=1)/11/g;
$idp  =~ s/10(?=1)/11/g;
# $lidp  =~ s/10(?=1)/11/g;
$gap  =~ s/10(?=1)/11/g;

# just an idea - how about smoothing those single column gaps away here?
$gap =~ s/01(?=0)/00/g;

$res = join('',@res);

output "res\t".$res."\n" ;
output "peak\t".$peak."\n";
output "idp\t".$idp."\n";
output "lidp\t".$lidp."\n";
output "gap\t".$gap."\n";

@peak = split(/ */, $peak);
@idp = split(/ */, $idp);
@lidp = split(/ */, $lidp);
@gap = split(/ */, $gap);

# 1. get peaks (2 in our case..)

my %status_ungapped;
# my $aln_len = $aln->length; 

foreach my $seq ($aln->each_seq) {

#    $status_gapped{ $seq->id } = '0'x$aln_len;

    my $seq_str = $seq->seq;
    $seq_str =~ s/-//g; # gapless..
    my $seq_len = length( $seq_str );
    my $id = $seq->id;
    $id =~ s/\s+$//;
    $status_ungapped{ $id } = '0'x$seq_len;
}

output "\nFind ultraconserved regions.\n\n";

my %motif_member_id;
my %motif_member_start; 
my %motif_member_end;

#my %motif_member_gapped_substr

my @used = (0) x scalar(@peak);
my $n_peaks = find_peak(@peak); # sets used..

my $found_motifs = scalar(keys %motif_member_id);
output "Found $found_motifs different motifs in $n_peaks ultraconserved regions.\n";

my @non_gapped_idp;

for (my $i = 0; $i < scalar(@idp); $i++) {
    # drop the peaks that were already blasted from idp 
    # they are probably clean (residual mop-up is left for the last stage)
    $non_gapped_idp[$i] = ( ($idp[$i]> 0) && ($gap[$i]==0) && ($used[$i]==0) ? 1 : 0);
}

output "\nFind ungapped high ID% regions, disregarding the previous peak motif parts.\n\n";

output "res\t".$res."\n";
output "used\t".(join("", @used))."\n";
output "ng_idp\t".(join("", @non_gapped_idp))."\n";

$n_peaks = find_peak(@non_gapped_idp);

my $current_motifs = scalar(keys %motif_member_id);
my $new_motifs = $current_motifs - $found_motifs;
$found_motifs = $current_motifs;

output "Found $new_motifs different motifs in $n_peaks different ungapped high ID% regions.\n";

# now, directly look at the co-aligned, but unhit seqs?

my @gapped_idp;
for (my $i = 0; $i < scalar(@idp); $i++) {
    $gapped_idp[$i] = ((($idp[$i]>0) && ($gap[$i]>0)) ? 1 : 0);
}

output "\nFind high ID% in gaps.\n\n";

output "res\t".$res."\n";
output "g_idp\t".(join("", @gapped_idp))."\n";
output "gap\t".(join("", @gap))."\n";

$n_peaks = find_peak(@gapped_idp);

$current_motifs = scalar(keys %motif_member_id);
$new_motifs = $current_motifs - $found_motifs;
$found_motifs = $current_motifs;

output "Found $new_motifs different motifs in $n_peaks different gapped high ID% regions.\n";

#my @non_gapped_lidp;

#for (my $i = 0; $i < scalar(@lidp); $i++) {
#    $non_gapped_lidp[$i] = (($lidp[$i]> 0) && ($gap[$i]==0) && ($used[$i]==0));        
#}

# $DEBUG && print "\nFind ungapped high ID% regions, disregarding the previous peak motif parts.\n\n";

output "\nLastly, find motifs in previously unhit parts.\n\n";

$DEBUG && output "res\t".$res."\n";
$DEBUG && output "peak\t".$peak."\n";
$DEBUG && output "ng_idp\t".(join("", @non_gapped_idp))."\n";
$DEBUG && output "g_idp\t".(join("", @gapped_idp))."\n";
$DEBUG && output "used\t".(join("", @used))."\n";

my $used= join("",@used);
$used =~ s/1/9/g;
$used =~ s/0/1/g;
$used =~ s/9/0/g;

my @untouched = split(/ */,$used);
find_peak(@untouched);

$current_motifs = scalar(keys %motif_member_id);
$new_motifs = $current_motifs - $found_motifs;

output "Found $new_motifs different motifs in $n_peaks different previously unhit regions.\n";

$DEBUG && output "\n\nStatus at end of loop:\n";
$DEBUG && output "res\t".$res."\n";

if($DEBUG) {
    foreach my $id ( keys %status_ungapped ) {
	output $id."\t".($status_ungapped{$id})."\n";
    }
}

# construct new indicator where ungapped is 

# so, did we cover the peak and id intervening sections properly?

# also, print UA*/UK* *S/*M ratios (based on number of actual reads or dominance only?)

my @motifs = (keys %motif_member_id);
#my @motif_seover = '0' x @motifs;
# withhold printing to allow pval-sort
my @motif_pval = 1 x @motifs;
my @motif_bulkout = "" x @motifs;

my $motif_nr = 0;
foreach my $seed_id ( @motifs ) {
    my %counted_isolate; 	    # new count for each motif
    
    my @motif_id = @{$motif_member_id{$seed_id}};
    my @motif_start = @{$motif_member_start{$seed_id}}; # motif positions converted to 0-base by save_motif
    my @motif_end = @{$motif_member_end{$seed_id}};

#    my $se_dominants_in_motif = 0;
#    my $co_dominants_in_motif = 0;
    my %isolatetype_dominants_in_motif;

    my $actual_members = scalar(@motif_id);

    # this is kind of silly, but we don't really want singletons, right?
    if( $actual_members < 2 ) {
	$DEVELOP_DEBUG && output "Not keeping motif $seed_id only $actual_members members.\n";
	next;
    }

    $motif_bulkout[$motif_nr] = "\nSeed $seed_id - $actual_members members\n";

    # print seed seq

    for (my $i = 0; $i < @motif_id; $i++) {
	
	my ($motif_id_seq) = ($aln->each_seq_with_id($motif_id[$i]));
	
	if (defined($motif_id_seq)) {
	    my $cand_str = $motif_id_seq->seq();
	    $cand_str =~ s/-//g;
	    my $cand_substr = substr( $cand_str, $motif_start[$i], $motif_end[$i]-$motif_start[$i]+1 );

	    $motif_bulkout[$motif_nr] .= sprintf "%10s\t%3d\t%3d\t%20s",$motif_id[$i],$motif_start[$i]+1, $motif_end[$i]+1, $cand_substr;
	} else {
	    $WARNING && output "WARNING: no seq by the id ".($motif_id[$i])."\n";
	    $motif_bulkout[$motif_nr] .= $motif_id[$i]."\t".($motif_start[$i]+1)."\t".($motif_end[$i]+1);
	}

	if(defined($fractcount{$motif_id[$i]}) && $fractcount{$motif_id[$i]} ne "") {
	    
	    my $fractcount = $fractcount{$motif_id[$i]};

	    if($count_one_dominant_per_isolate == 1) {	

		# isolate name case dependent - feature?
		while ($fractcount =~ m/([a-zA-Z]+)(\d+)/g ) {

		    my $isolate_type = $1;

		    if ($SE_UKS_COMPAT) {
			# old, presumably redundant check for error in naming convention
			if( $isolatetype eq "cp" ) {
			    $isolatetype = "co";
			}

			# group uks & se for final tally
			if ( $isolatetype eq "uks" ) {
			    $isolatetype = "se";
			}
		    }

		    my $current_count = $isolatetype.$2;

		    if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {
		    } else {	    

			$counted_isolate{$current_count} = 1;

			$isolatetype_dominants_in_motif{$isolatetype}++;
		    }
		}

	    } else {
		
		while ($fractcount =~ m/([a-zA-Z]+)\d+/g ) {
		    
		    my $isolatetype = $1;
		    
		    if ($SE_UKS_COMPAT) {

			# old, presumably redundant check for error in naming convention	       
			if( $isolatetype eq "cp" ) {
			    $isolatetype = "co";
			}

			# group uks & se for final tally
			if ( $isolatetype eq "uks" ) {
			    $isolatetype = "se";
			}
		    }

		    $isolatetype_dominants_in_motif{$isolatetype}++;
		}

	    }

	    $motif_bulkout[$motif_nr] .= "\t".$fractcount;

	}
	$motif_bulkout[$motif_nr] .= "\n";
    }
    
    my (@isolatetypes) = keys %total_dominants_in_motif;

    if( @isolatetypes >2) {
	$WARNING && output "WARNING: more than two (unaggregated) isolatetypes found in sequence set. Initially comparing the two first, ", join(", ", @isolatetypes[0..1]),". Please set up separate dominance tables and rerun to score other combinations.\n";
    }

    ### EDIT POINT: some remain for OPI-model above, and far below for dominance tab parser, Ok.
    ### as well as the rather important SE dominants: \d+ CO dominants: \d+ line, printed from this main, far above, at dominance tab parser call
    # implement direct call score_using_R?

    my $se_type = $isolatetypes[0];
    my $co_type = $isolatetypes[1];
    
    if($isolatetype_dominants_in_motif{$se_type} + $isolatetype_dominants_in_motif{$co_type} > 0 ) {
	
	$DEBUG && output "DEBUG: Calling dhyper for seed $seed_id ($actual_members actual members) with M=".($isolatetype_dominants{$se_type}).", N=".($isolatetype_dominants{$co_type}).", m=".($isolatetype_dominants_in_motif{$se_type}).", n=".($isolatetype_dominants_in_motif{$co_type})."\n";
	my ($pval, $se_over) = dhyper($isolatetype_dominants{$se_type}, $isolatetype_dominants{$co_type}, $isolatetype_dominants_in_motif{$se_type}, $isolatetype_dominants_in_motif{$co_type});
	$motif_pval[$motif_nr] = $pval;
#	$motif_seover[$motif_nr] = $se_over;

	$motif_bulkout[$motif_nr] .= sprintf "%s %3d %s %3d %sfraction %.3f p %.4f %s\n", uc($se_type), $isolatetype_dominants_in_motif{$se_type}, uc($co_type), $isolatetype_dominants_in_motif{$co_type}, uc($se_type), ($isolatetype_dominants_in_motif{$se_type} / ( $isolatetype_dominants_in_motif{$se_type}+$isolatetype_dominants_in_motif{$co_type})), $pval, $se_over;
#	exec("echo |R")
    }
    $motif_nr++;
}

print @section_log;

print "\nMotifs\n";
# sort
print "Found $motif_nr motifs.\n";

my @order;
#if ($fractioncount_file ne "") {

@order = sort {$motif_pval[$a] <=> $motif_pval[$b] } (0..($motif_nr-1));   

#} else {
#    @order = 0..($motif_nr-1);
#}

foreach my $i ( @order ) {
    print $motif_bulkout[$i];
}

$DEBUG && print "\nDone. ;-)\n";

sub subseq_from_columns {

    my $seed_seq = shift;
    my $peak_start = shift;
    my $peak_end = shift;

    my $peak_seq_loc_start = $seed_seq->location_from_column($peak_start);
    if (!defined($peak_seq_loc_start) || $peak_seq_loc_start->location_type() eq 'IN-BETWEEN') {
	$DEBUG && output "Sequence ".($seed_seq->id)." has no perfect match (at start)..\n";
	return ("", 0,0,0,0);
    }
    my $peak_seq_start = $peak_seq_loc_start->start; # 1-based inclusive
    
    my $peak_seq_loc_end = $seed_seq->location_from_column($peak_end);
    if (!defined($peak_seq_loc_end) || $peak_seq_loc_end->location_type() eq 'IN-BETWEEN') {
	$DEBUG && output "Sequence ".($seed_seq->id)." has no perfect match (at end)..\n";
	return ("", 0,0,0,0);
    }
    my $peak_seq_end = $peak_seq_loc_end->end; # 1-based inclusive
    
    my $rep_cand_str = $seed_seq->seq();
    $rep_cand_str =~ s/-//g;
    my $rep_cand_len = length($rep_cand_str);

    my $peak_len = $peak_seq_end-$peak_seq_start+1;

    return ($rep_cand_str, $rep_cand_len, $peak_seq_start, $peak_seq_end, $peak_len);
}

sub find_peak {

    my @indicator = @_;
    my $n_peaks = 0;

    # my $res = shift;
    
    my $peak_start = -1;
    my $peak_end = -1;

    for (my $i = $aln_start; $i <= $aln_end; $i++) {

	if ($peak_start == -1) {
	    $DEBUG && output "Examining pos $i; not at peak at the moment.\n";
	    if($indicator[$i-1] == 1) {
		$peak_start = $i;
		$peak_end = $i;
	    }
	} else {
	    
	    if(($indicator[$i-1] == 0) || ($i == $aln_end)) {
		# leaving peak

		if( $peak_start == $peak_end ) {
		    $DEBUG && output "Ignoring single peak residue at pos $peak_start.\n";

		    $peak_start= -1;
		    $peak_end= -1;
		    next;
		}

#		output  "Found indicator peak between $peak_start and $peak_end.\n";

		output "Indicator peak: $peak_start - $peak_end.\n";
		$n_peaks++;

		# some of this should be generalised to fit for check_remaining_hits_in_region as well..

		# 1. pick a peak-representative read
		
		my $seed_seq;
		my $found_representative = 0;
		
		my $peak_string = substr($res, $peak_start-1, $peak_end-$peak_start+1);
		
		my $lext = 0;
		my $rext = 0;

		$DEBUG && output "Looking for representative for $peak_string.. ($peak_start - $peak_end)\n";
		for(my $j = 1; $j <= $aln_height; $j++) {

		    $seed_seq = $aln->get_seq_by_pos($j);
		    my ($rep_cand_str, $rep_cand_len, $peak_seq_start, $peak_seq_end, $peak_len) = subseq_from_columns($seed_seq, $peak_start, $peak_end);
		    if ($rep_cand_str eq "") {
			next; # failed; only gaps around here?
		    }
		    $rep_cand_substr = substr($rep_cand_str, $peak_seq_start-1, $peak_len); # substr takes 0-based coord, peak coord is 1-based
		    $rep_cand_substr_len = length($rep_cand_substr);

		    # silly cautious: no extension possible if entire sequence of tentative seed is shorter than min_blastp_len
		    if ($rep_cand_len < $min_blastp_len) {
			$WARNING && output "Candidate contig ".($seed_seq->id)." is shorter $rep_cand_str than min blastp query ($min_blastp_len).\n";
			next;
		    }

		    $DEBUG && output "Sequence ".($seed_seq->id)." $peak_seq_start - $peak_seq_end : $rep_cand_substr";

		    if( $rep_cand_substr eq $peak_string ) {
			$DEBUG && output "..match!\n";

			# Extend hit if too short for blast
			if($rep_cand_substr_len < $min_blastp_len) {			    
			    
			    if ($status_aware) {			    
				($peak_len, $peak_seq_start, $peak_seq_end, $lext, $rext) = extend_hit_region_status_aware($peak_len, $peak_seq_start, $peak_seq_end, $rep_cand_str, $rep_cand_len, $seed_seq->id);
			    } else {
				($peak_len, $peak_seq_start, $peak_seq_end, $lext, $rext) = extend_hit_region($peak_len, $peak_seq_start, $peak_seq_end, $rep_cand_str, $rep_cand_len);
			    }

			    $rep_cand_substr = substr($rep_cand_str, $peak_seq_start-1, $peak_seq_end-$peak_seq_start+1);
			    $rep_cand_substr_len = length( $rep_cand_substr );
#			    $peak_len = $peak_seq_end-$peak_seq_start+1; # duh?

			    $WARNING && output "(Peak ($peak_string) too short for blastp; using ".($seed_seq->id).
			    " $peak_seq_start - $peak_seq_end ($rep_cand_substr) as a representative.\n";
			}

			# check status! dont bother calling with the unextended motif - the rules for overlaps are not very compatible with shorts..
			my $tried = check_domain_status($seed_seq->id, $peak_seq_start, $peak_seq_end); # call with 1-based incl.

			if($tried > 0) {
			    $WARNING && output "Candidate contig ".( $seed_seq->id )." $peak_seq_start - $peak_seq_end has status $tried. Not a rep.\n";
			    next;
			}

			if ($rep_cand_substr_len < $min_blastp_len) {
			    $WARNING && output "Candidate ".( $seed_seq->id)." is shorter $rep_cand_substr than min blastp query ($min_blastp_len) after extension.\n";
			    next;
			}

			$found_representative = 1;

			# mark query as tested -- really not needed.. Should hopefully be hit by its own representative.. =)
			my $save = markup_domain($seed_seq->id, $peak_seq_start, $peak_seq_end, 1);

			# also mark this region in consensus as "used"; so far only used to remove "peaks" from "idp"-regions
			@used[($peak_start-1)..($peak_end-1)] = (1) x ($peak_end - $peak_start+1);
			
			# construct representative for blasting
			my $old_id = $seed_seq->id;
			$old_id =~ s/\s+$//;
			$seed_seq = Bio::Seq->new( '-id' => $old_id."_".$peak_seq_start."_".$peak_seq_end, '-seq' => $rep_cand_substr);
			if($save == 1) {
			    save_motif($seed_seq->id, $old_id, $peak_seq_start, $peak_seq_end);
			    $DEBUG && output "Saving motif ".($seed_seq->id).", $old_id, $peak_seq_start, $peak_seq_end [find_peak]\n";
			}

			last;
		    } else {
			$DEBUG && output "..mismatch!\n";
		    }
		}
		
		if ($found_representative == 0) {
#		    $WARNING && output "Could not find exact representative for peak at columns $peak_start - $peak_end ($peak_string). Proceeding to next motif instead.\n";
		    $WARNING && output "Could not find exact representative for peak at columns $peak_start - $peak_end ($peak_string). Proceeding to by row motif search instead.\n";
		    # next;
		    
		} else {		    
		    # 2. run blastpgp
		    run_blastpgp($seed_seq, $peak_start-$lext, $peak_end+$rext);
		}

		$DEBUG && output "\nCheck remaining hits in columns $peak_start - $peak_end..\n";
		check_remaining_hits_in_region($peak_start, $peak_end); # order of this can be rather important, really. But, I recall the first versions did not use check_remaining. Given good hit markup, this should be increasing some sort of sensitivity, but slow
		# And coords - extend or not? FOLLOW-UP 

# then set used: c'mon we've been through every seq 
		@used[($peak_start-1)..($peak_end-1)] = (1) x ($peak_end - $peak_start+1);

		$DEBUG && output "\nLeft checking remaining hits $peak_start - $peak_end..\n";

		$peak_start= -1;
		$peak_end= -1;

		# then, immediately check any remaining non-hits as queries, until all have been tried
		
	    } elsif($indicator[$i-1] == 1) {
		$peak_end = $i;
	    }
	}
    }
    return $n_peaks;
}

sub extend_hit_region {
    my $peak_len = shift;
    my $peak_seq_start = shift;
    my $peak_seq_end = shift;
    my $rep_cand_str = shift;
    my $rep_cand_len = shift;

    my $in_start = $peak_seq_start; 
    my $in_end = $peak_seq_end; 

    my $needed_extension = $min_blastp_len - $peak_len;
    
    my $needed_extension_first_half = POSIX::ceil( $needed_extension / 2 );
    my $needed_extension_last_half = POSIX::floor( $needed_extension / 2 );

    if( $peak_seq_start-1 - $needed_extension_first_half < 0 ) {
	$needed_extension_last_half += $needed_extension_first_half - ($peak_seq_start - 1);
	$peak_seq_start = 1;
    } else {
	$peak_seq_start -= $needed_extension_first_half;
    }

    if( $peak_seq_end + $needed_extension_last_half > $rep_cand_len ) {
	$peak_seq_start -= $needed_extension_last_half - ($peak_seq_end - $rep_cand_len);
	($peak_seq_start < 1) && ($peak_seq_start = 1);
	$peak_seq_end = $rep_cand_len;
    } else {
	$peak_seq_end += $needed_extension_last_half;
    }

    # update candidate substrin
    $rep_cand_substr = substr($rep_cand_str, $peak_seq_start-1, $peak_seq_end-$peak_seq_start+1);
    $peak_len = $peak_seq_end-$peak_seq_start+1;

    my $left_ext = $in_start - $peak_seq_start;
    my $right_ext = $peak_seq_end - $in_end;

    return ($peak_len,$peak_seq_start, $peak_seq_end, $left_ext, $right_ext);
}

sub extend_hit_region_status_aware {
    my $peak_len = shift;
    my $peak_seq_start = shift;
    my $peak_seq_end = shift;
    my $rep_cand_str = shift;
    my $rep_cand_len = shift;

    my $id = shift;

    # save original coordinates
    my $in_start = $peak_seq_start; 
    my $in_end = $peak_seq_end; 

    my $left_ext = 0;
    my $right_ext = 0;

    # attempt the usual extension
    ($peak_len, $peak_seq_start, $peak_seq_end, $left_ext, $right_ext) = extend_hit_region($peak_len, $peak_seq_start, $peak_seq_end, $rep_cand_str, $rep_cand_len);    

    # update candidate substring
    $rep_cand_substr = substr($rep_cand_str, $peak_seq_start-1, $peak_seq_end-$peak_seq_start+1);
    
    $DEVELOP_DEBUG && output "status_aware_extension: original extension suggestion of ($in_start,$in_end) is ($peak_seq_start,$peak_seq_end) ($rep_cand_substr)\n";
    my $status = substr($status_ungapped{$id}, $peak_seq_start, $peak_len); #oops - 0based?
    $DEVELOP_DEBUG && output "status original extension $status, ";
    
    if ($status ne ('0' x $peak_len )) {

	# give us a shove in the direction away from already hit/used parts, if thats simple to determine.
#	my (@initial_set_status) = ($status =~ m/[12]{1}/g);
	
	if ($status=~m/^([12]+)0+/) {
	    # attempt to extend further to the right..
	    my $left_olap = length($1);

	    $DEVELOP_DEBUG && output "so shove to the right.\n";
	    if($peak_seq_end + $left_olap <= $aln_end) {
		$peak_seq_end += $left_olap;
		# also take in the the left olap..
		$peak_seq_start += $left_olap;
	    } elsif( $peak_seq_end < $aln_end ) { 

		my $diff = $aln_end - $peak_seq_end;
		
		$peak_seq_end = $aln_end;
		$peak_seq_start += $diff;
	    }
	} elsif ($status =~ /0+([12]+)$/) {
	    # attempt to shift extension to the left
	    my $right_olap = length($1);

	    $DEVELOP_DEBUG && output "so shove to the left.\n";	   
	    if($peak_seq_start - $right_olap > $aln_start) {
		$peak_seq_end -= $right_olap;
		$peak_seq_start -= $right_olap;
	       

	    } elsif( $peak_seq_start > $aln_start ) {
		my $diff = $peak_seq_start - $aln_start;
		$peak_seq_start = $aln_start;
		$peak_seq_end -= $diff;
	    }	    
	}
    
        $rep_cand_substr = substr($rep_cand_str, $peak_seq_start-1, $peak_seq_end-$peak_seq_start+1);
	$peak_len = $peak_seq_end-$peak_seq_start+1;	
    }
    
    $status = substr($status_ungapped{$id}, $peak_seq_start, $peak_len); #oops - 0based?

    $DEVELOP_DEBUG && output "Currrent extension suggestion is ($peak_seq_start,$peak_seq_end) ($rep_cand_substr), status $status.\n";

    $left_ext = $in_start - $peak_seq_start;
    $right_ext = $peak_seq_end - $in_end;
    
    return ($peak_len,$peak_seq_start, $peak_seq_end, $left_ext, $right_ext);
}

sub extend_library_region {

    my $peak_seq_start = shift;
    my $peak_seq_end = shift;    

    my $peak_len = $peak_seq_end - $peak_seq_start + 1;
    my $needed_extension = 3 * $peak_len ;
#    my $rep_cand_str = shift;
#    my $rep_cand_len = shift;
    
    my $needed_extension_first_half = POSIX::ceil( $needed_extension / 2 );
    my $needed_extension_last_half = POSIX::floor( $needed_extension / 2 );

    # EXTEND_DEBUG
#    print STDERR "Needed_ext: $needed_extension\n";    
#    print STDERR "first half: $needed_extension_first_half second half: $needed_extension_last_half\n";
        	    
    if( $peak_seq_start-1 - $needed_extension_first_half < 0 ) {
	$needed_extension_last_half += $needed_extension_first_half - ($peak_seq_start - 1);
	$peak_seq_start = 1;
#	print STDERR "C1: needed_ext_last_half now $needed_extension_last_half\n";
    } else {
	$peak_seq_start -= $needed_extension_first_half;
#	print STDERR "C1e: peak_seq_start: $peak_seq_start\n";
    }
    
    if( $peak_seq_end + $needed_extension_last_half > $aln_end ) {
	$peak_seq_start -= $needed_extension_last_half - ($peak_seq_end - $aln_end);
	($peak_seq_start) < 1 && ($peak_seq_start = 1);
	$peak_seq_end = $aln_end;
#	print STDERR "C2: peak_seq_start: $peak_seq_start, end at end $aln_end\n";
    } else {
	$peak_seq_end += $needed_extension_last_half;
#	print STDERR "C2e: peak_seq_end: $peak_seq_end\n";
    }
    
    # update candidate substrin
#    $rep_cand_substr = substr($rep_cand_str, $peak_seq_start-1, $peak_seq_end-$peak_seq_start+1);
#    $peak_len = $peak_seq_end-$peak_seq_start+1;

    return ($peak_seq_start,$peak_seq_end);
}

sub check_remaining_hits_in_region {

    my $region_start = shift;
    my $region_end = shift;

  ROW:    for(my $j = 1; $j <= $aln_height; $j++) {

      my $seed_seq = $aln->get_seq_by_pos($j);
      
      my ($seed_seq_str, $seed_seq_len, $seed_seq_start, $seed_seq_end, $peak_len) = subseq_from_columns($seed_seq, $region_start, $region_end);
      if ($seed_seq_str eq "") {
	  next; # failed, presumably due to gap at start/end of interval - leave for final mop-up?
      }
      my $seed_seq_substr = substr($seed_seq_str, $seed_seq_start-1, $seed_seq_end-$seed_seq_start+1);
      my $seed_seq_substr_len = length($seed_seq_substr);
      # check status!

      # Extend hit if too short for blast
	
      # silly cautious: no extension possible if entire sequence of tentative seed is shorter than min_blastp_len
      if ($seed_seq_len < $min_blastp_len) {
	  $WARNING && output "Candidate ".( $seed_seq->id ). " is shorter $seed_seq_str than min blastp query ($min_blastp_len).\n";
	  next; 
      }

      my $rext = 0;
      my $lext = 0;

      if($seed_seq_substr_len < $min_blastp_len) { # difficult: use seed-seqlen or peak-len?
	  $DEBUG && output "Attempting to extend...\n";
	  ($peak_len, $seed_seq_start, $seed_seq_end, $lext, $rext) = extend_hit_region($peak_len, $seed_seq_start, $seed_seq_end, $seed_seq_str, $seed_seq_len);

	  $seed_seq_substr = substr($seed_seq_str, $seed_seq_start-1, $seed_seq_end-$seed_seq_start+1);
	  $seed_seq_substr_len = length($seed_seq_substr);

	  $WARNING && output "Peak too short (now $peak_len) for blastp; using ".( $seed_seq->id).
	  " $seed_seq_start - $seed_seq_end ($seed_seq_substr) as a representative.\n";
      }

      my $tried = check_domain_status($seed_seq->id, $seed_seq_start, $seed_seq_end);
	
      if($tried > 0) {
	  $WARNING && output "Candidate contig ".( $seed_seq->id )." $seed_seq_start - $seed_seq_end has status $tried. No need for additional search!\n";
	  next ROW; # get next row instead
      }

      $DEBUG && output "\nCheck remaining hits says sequence ".($seed_seq->id)." $seed_seq_start - $seed_seq_end : $seed_seq_substr was untouched.\n";

      # mark query as tested 
      my $save = markup_domain($seed_seq->id, $seed_seq_start, $seed_seq_end, 1);
      # no need to save status -- really not needed.. Should hopefully be hit by its own representative.. =)

      my $old_id = $seed_seq->id;
      $old_id =~ s/\s+$//;
      $seed_seq = Bio::Seq->new( '-id' => $old_id."_".$seed_seq_start."_".$seed_seq_end, '-seq' => $seed_seq_substr);

      if($save == 1) {
	  save_motif($seed_seq->id, $old_id, $seed_seq_start, $seed_seq_end);
	  $DEBUG && output "Saving motif ".($seed_seq->id).", $old_id, $seed_seq_start, $seed_seq_end [from check remaining hits]\n";
      }

      # ehhhrm, we don't want to see any singletons, do we?!

      # nb - the extended sequence coordinates can extend outside of the region
      # AND Im not convinced they are obi-one proof - please recheck!
      run_blastpgp($seed_seq, $region_start-$lext, $region_end+$rext); # sending entire region coords for the mop-up phase..


    }
}

sub run_blastpgp {

    my $seq = shift;
    # my $database = shift; # first, make sure there is a database.. 
    # not currently adaptive. using global.
    my $seed_region_start =shift; # gapped coords
    my $seed_region_end = shift;  # gapped coords
    my $seq_name = $seq->id;
    
    # ungapped offsets for alignment subsection when using strict sectioning
    my $seq_relative_region_start = 0;
#    my $seq_relative_region_end = 0;

    # get sub-part of alignment, and press new database
    # need to keep track of start (end?) of rows to map coords back to the original aln for markup
    #   -- 1 coord for each start: the blastentries are gapless, aln is gapped    

   #  1: print subpart to see that we get it right.
    
    my $local_db_file = $seq->id.".db.fasta"; # ehhrr and else!!
    my $outfile_basename = $local_db_file."_searched_by_".$seq_name;

    # use one total file iff no strict
    if( $strict_section == 0 ) {

	$local_db_file = $db_file; # .".db.fasta"; # ehhrr and else!!
	$outfile_basename = $db_file."_searched_by_".$seq_name;

    } else {

	$local_db_file = $seq->id.".db.fasta"; # ehhrr and else!!
	$outfile_basename = $local_db_file."_searched_by_".$seq_name;

#	open CLEAR, ">".$local_db_file.".sn";
#	close CLEAR;

#	open CLEAR, ">".$local_db_file.".pn";
#	close CLEAR;	
    
    }

    # for saving strict offset for each aln member in this section
    my %seq_rel_region_start;
#    my %seq_rel_region_end;

    if ( $strict_section == 1) {
	
	my $db_temp_seqio = Bio::SeqIO->new('-file'=>">$local_db_file", '-format'=>'Fasta');

	$DEVELOP_DEBUG && output "Strict section from $seed_region_start to $seed_region_end for seed ".($seq->id).".\n";
#	$DEVELOP_DEBUG && print STDERR  "Strict section from $seed_region_start to $seed_region_end for seed ",$seq->id,".\n";

	($seed_region_start,$seed_region_end) = extend_library_region($seed_region_start, $seed_region_end);

	$DEVELOP_DEBUG && output "Strict extends to $seed_region_start - $seed_region_end.\n";
#	$DEVELOP_DEBUG && print STDERR "Strict extends to $seed_region_start - $seed_region_end.\n";
	
ALNSEQ:	foreach my $alnseq ($aln->each_seq) {
	    # pick coords before removing gaps
	    my $seq_str = substr($alnseq->seq(),$seed_region_start-1,$seed_region_end-$seed_region_start+1); # aln coords are 1-based inclusive
	    my $seq_id = $alnseq->id;
	    $seq_id =~ s/\s+$//;
	    
	    $DEVELOP_DEBUG && output "Strict: $seq_str\t$seq_id.\n";

	    $seq_str =~ s/-//g;
	    if( length($seq_str) == 0  ) {
		$DEVELOP_DEBUG && output "Strict: oopsidaisy - no aa in $seq_id in aln - not adding this to db.\n";
		next ALNSEQ;
	    } elsif( length($seq_str) < $min_blastp_len) {
		$DEVELOP_DEBUG && output "Strict: oopsidaisy - fewer than threshold $min_blastp_len aa in $seq_id in aln; not adding to db.\n";
		next ALNSEQ;
	    }

	    # save per-seq ungapped coords to be able to accurately mark the returned hits..
	    my $alnseq_start_loc = $alnseq->location_from_column($seed_region_start);
	    if (!defined($alnseq_start_loc)) {
#		output "Abunai! No section start defined for $seq_id.  Assuming that all undef starts are prior to actual sequence, and setting start to 1.\n";
		$seq_relative_region_start = 1;
	    } elsif ($alnseq_start_loc->location_type() ne 'exact') {
		$seq_relative_region_start = $alnseq_start_loc->max_start;
	    } else {
		$seq_relative_region_start = $alnseq_start_loc->max_start;
		$DEVELOP_DEBUG && output "Start min $seq_relative_region_start for ".($alnseq_start_loc->location_type())."\n";
	    }

#	    my $seq_relative_region_end = $alnseq->location_from_column($seed_region_end)->max_end;
	    
#  save as features..
#	    $seq->add_SeqFeature(new Bio::SeqFeature::Generic ( -start => $seq_relative_region_start, 
#								-end => $seq_relative_region_end, 
#								-strand => 1, 
#								-primary => 'seed_loc',
#								-tag => {
#								    aln_start => $aln_rel_start,
#								    aln_end => $aln_rel_end
#								    }
#								) );
	    # or just drop them in a hash on seq_id? less elegant, but efficient.
	    $seq_rel_region_start{$seq_id} = $seq_relative_region_start;
	    
#	    $seq_rel_region_end{$seq_id} = $seq_relative_region_end;
	    # should store those in a local hash, and send to the resul

	    $outseq = Bio::Seq->new( '-id' => $seq_id, '-seq' => $seq_str );
	    $db_temp_seqio->write_seq($outseq);
	}
	
	system("formatdb -i $local_db_file -p T");
	$DEVELOP_DEBUG && output "Strict section: pressed region for seed ".($seq->id).".\n";
    } else {
	# use normal db..
    }

    # write the sequence file... filename .ceq
    my $ceq_file = $outfile_basename.".ceq";

    open CEQ, ">$ceq_file";
    print CEQ ">$seq_name\n";
    print CEQ $seq->seq(),"\n";
    close CEQ;

    # sn and pn are global

    # add ceq-file name to list (kind of slow to use append; the alternative is to keep an open filehandle from caller)
    open SNLIST, ">>$db_file.sn";
    print SNLIST $ceq_file."\n";
    close SNLIST;

    my $chk_file = $outfile_basename.".blpgp.chk"; 

    # add the to be produced chk-file name to list (kind of slow to use append; the alternative is to keep an open filehandle from caller and pass it, or use a global glob)
    open PNLIST, ">>$db_file.pn";
    print PNLIST $chk_file, "\n";
    close PNLIST;

    @params = ('database' => $local_db_file ,'outfile' => $outfile_basename.".blpgp",
	       'e' => $blastpgp_e, 'h'=>$blastpgp_h, 'f'=>$blastpgp_f, 'F'=>$blastpgp_F, 'b'=>$blastpgp_b, 'v'=>$blastpgp_v,
	       # testade e=1000 med dåliga resultat..
#	       'm'=>6,
	       'C'=>$chk_file,
	       '_READMETHOD'=> 'BPlite'
#	       '_READMETHOD'=> 'Blast' #seems fairly good, but requires quite a few changes to the iteration below..
	       );
# m=6 is a flag to get aln working according to BPlite perldoc
    $DEBUG && output "Blastpgp of $seq_name (".($seq->seq).") vs $local_db_file..\n";

    my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
    $factory->j(10);
    my $report = $factory->blastpgp($seq);

    if( $report->{'REPORT_DONE'} == 1 ) {
	$WARNING && output "OOPS: Report reported done before calling number_of_iterations.\n";
	$WARNING && output "This probably means there were no hits below e-val treshold. Dropping current query.\n";
	return;
	# note that the seed has already been marked status 1..
    }

    my $total_iterations = $report->number_of_iterations();
    my $last_iteration = $report->round($total_iterations);
   
#    my $result = $report->next_result();  #Blastparsetest#
    
#    my $last_iteration = $result->iteration();  #Blastparsetest# 

    # DEBUGish alignment dump -- only works for m==6; HSPs do NOT work properly for m==6 at the moment
# $factory->m(6);
# $report = $factory->blastpgp($seq);
# $total_iterations = $report->number_of_iterations;
# $last_iteration = $report->round($total_iterations);
#    my $align = $last_iteration->Align;
#    if ($align->no_sequences < 2) {
	# drop this alignment/motif
#	$WARNING && print "Number of hits < 2, so dropping current query.\n";
#	return;
#    }
    
#    my $out_aln = Bio::AlignIO->new('-format' => 'clustalw', '-file' => ">".$outfile_basename.".psi.aln");    
#    $out_aln->write_aln($align); # names broken!? Why?

#    my $result = 

#    hits_below_threshold #blastparsetest#

    my $hsp_pos_cutoff = $last_iteration->qlength - $hsp_neg_cutoff;
    
    my @markupqueue;
    while (my $sbjct = $last_iteration->nextSbjct) {
	
	$DEBUG && output "Iterated blast hit ", $sbjct->name,".\n";
	
	while (my $hsp = $sbjct->nextHSP) {
	    if($hsp->positive() >= $hsp_pos_cutoff) {
		$DEBUG && output sprintf "Motif member %10s  %3d - %3d : %s\n", $hsp->hit->seq_id, $hsp->hit->start, $hsp->hit->end, $hsp->sbjctSeq;
		my $hsp_hit_start = $hsp->hit->start;
		my $hsp_hit_end = $hsp->hit->end;

		if($strict_section == 1) {
		    my $seq_id = $hsp->hit->seq_id;
		    $seq_id =~ s/\s+$//g;
# $seq_id_seed =~ m/[^_]+)/ );

		    my $offs = $seq_rel_region_start{$seq_id} -1; # aln start (stored in seq_rel_reg_start) is 1-based, as is hsp hit start
		    $DEVELOP_DEBUG && output "Strict section: Adjusting hit offset for $seq_id by $offs.\n";
#		    $seq_rel_region_end{$seq_id} = $seq_relative_region_end;
		    
		    $hsp_hit_start += $offs ;
		    $hsp_hit_end += $offs ;
		}

		push @markupqueue, [$hsp->hit->seq_id, $hsp_hit_start, $hsp_hit_end, 2];
                #
		# DEVOLP_DEBUG
		my $debug_motif_seq_id = $hsp->hit->seq_id;
		$debug_motif_seq_id =~ s/\s+$//;
	       
		$DEVELOP_DEBUG && output "Retrieving sequence $debug_motif_seq_id.\n";
		my ($debug_motif_seq) = ($aln->each_seq_with_id($debug_motif_seq_id));
		my $debug_motif = $debug_motif_seq->seq();
		$debug_motif =~ s/-//g;
		my $debug_motif_substr = substr( $debug_motif, $hsp_hit_start-1, $hsp_hit_end - $hsp_hit_start + 1 );

		$DEVELOP_DEBUG && output "Strict section: offset results in motif str ".$debug_motif_substr.".\n";
	    } else {
		$DEBUG && output "Dropping HSP since positives ".($hsp->positive)." less than $hsp_pos_cutoff (".($hsp->sbjctSeq).") .\n";
		# note that this is very powerless, since later motif identification will anyway use all hsps in for the matrix from the .chk-file!?! Better be open about it?
		# still useful if all but one have too many mismatches.. 
	    }
	}
    }
    if (@markupqueue > 0) { # the seed is already saved, so the first hit is the second member
	# mark subject hit!
	while (my $markme = shift @markupqueue) {
	    my $save = markup_domain(@$markme);
	    if($save == 1) {
		# RECHECK - is this correct?
		save_motif($seq_name, @$markme);
		$DEBUG && output "Saving motif $seq_name ".(join(" ",@$markme))." from run_blastpgp for seed ".$$markme[0]."\n";
	    }
	}
    } else {
	# drop this alignment/motif
	$WARNING && output "Number of subjects hit < 1, so dropping current query.\n";
	# note that the seed has already been marked status 1.. so don't print it..	
    }
    if($strict_section == 1) {
	system("rm $local_db_file.phr $local_db_file.pin $local_db_file.psq") unless $DEBUG;
    }
}

sub save_motif {

    my $seed_id = shift;
    my $id = shift;
    my $start = (shift) - 1; #convert to internal 0-based desc.
    my $end = (shift) - 1;
    my $type = shift; #unused here, only in markup domain

    $seed_id =~ s/\s+$//;
    $id =~ s/\s+$//;

    push @{$motif_member_id{$seed_id}}, $id;
    push @{$motif_member_start{$seed_id}}, $start;
    push @{$motif_member_end{$seed_id}}, $end;

    return;
}

sub markup_domain {
    # basically a hash of strings, where each string holds 
    # 0: untested, unhit
    # 1: tested
    # 2: hit by other element

    my $id = shift;
    my $start = (shift) - 1; #convert to internal 0-based desc.
    my $end = (shift) - 1;
    my $type = shift;

    my $hit_len = $end - $start + 1;

    my $tolerated_overlap_with_previous_hits = $tolerated_motif_overlap;

    $id =~ s/\s+$//;

    my $prev = substr($status_ungapped{$id}, $start, $hit_len);
    if ($prev ne ('0' x $hit_len )) {
	my (@prev_hits) = ($prev =~ m/[12]{1}/g); 
	if (@prev_hits > $tolerated_overlap_with_previous_hits) {
	    $WARNING && output "WARNING: Extensive double handling of $id $start - $end: was $prev is ".($type x $hit_len).". Ignoring match.\n";
	    return 0;
	}
    }
    # actuall replace

    $DEBUG && output "setting status for $id $start - $end to ".($type x $hit_len).".\n";
    substr($status_ungapped{$id}, $start, $hit_len, "$type" x $hit_len);
    return 1;
}

sub check_domain_status {
    # basically a hash of strings, where each string holds 
    # 0: untested, unhit
    # 1: tested
    # 2: hit by other element

    my $id = shift;
    my $start = (shift) - 1; #convert to internal 0-based desc.
    my $end = (shift) - 1;

    my $len = $end - $start + 1;

    $id =~ s/\s+$//;
    
    my $status = substr($status_ungapped{$id}, $start, $len);

    if ($status ne ('0' x $len )) {
	my (@set_status) = ($status =~ m/[12]{1}/g);
	
	if (@set_status > $tolerated_motif_overlap) {
	    $WARNING && output "WARNING: Previously touched $id $start - $end: $status. Ignoring potential seed.\n";
	    return(1);
	}

	if( $len < $tolerated_motif_overlap && @set_status > 0) {
	    # if overlap too short tolerance becomes stupid.. at tolerated motif overlap of 3, 1 previous used aa feels enough.
	    # at four, 2 should perhaps be enough
#	    if( @set_status == $len  ) {
		
		# the entire, short, motif is set - ignore.
	    $WARNING && output "WARNING: Previously touched $id $start - $end: $status. Ignoring potential seed.\n";
	    return(1);
#	    }    # 11100 1110 11100
	} elsif ( $len < 2 * $tolerated_motif_overlap && @set_status > 1 ) {
	    $WARNING && output "WARNING: Previously touched $id $start - $end: $status. Ignoring potential seed.\n";
	    return(1);
	}

	$WARNING && output "WARNING: Previously touched $id $start - $end: $status. Using anyway.\n";
	return(0);
    } else {
	$DEBUG && output "DEBUG Checking status for $id $start - $end: $status. Untouched.\n";
	return(0);
    }
}

#sub get_unchecked_parts {
    # will we really get any? a prining diagnostic is probably more useful
#    my $id = shift;
# my @status = split (/ */, $status_ungapped{$id});
    #my @seq = $aln->get_seq_by_id($id); # these alignment names *shoud* be unique
    #if (@seq > 1) {
#	$WARNING && print "Warning! Got ",scalar(@seq), " sequences with id $id from alignment when trying to check for unused parts.\n";
#    }
 #   my $seqstr = $seq[0]->seq();
    # check status for any longer than cutoff (6) consecutive unmapped regions    
#}

sub read_dominance_table {

    my $dt = shift; 
    
    while (my $row = <$dt>) { # $$dt ?
	my ($contig,$dominance_str) = split(/\s+/,$row);

	$DEBUG && output "$contig has dominance $dominance_str.\n";

	$fractcount{$contig} = $dominance_str; 
    }
    
    # ehh, no OPI? feature??
    foreach my $contig (keys %fractcount) {

	my $fractcount = $fractcount{$contig};

	while ($fractcount =~ m/([a-zA-Z]+)\d+/g ) {
	    
	    my $isolatetype = $1;
	
	    if ($SE_UKS_COMPAT) {
	    
		# old, presumably redundant check for error in naming convention	       
		if( $isolatetype eq "cp" ) {
		    $isolatetype = "co";
		}
		
		# group uks & se for final tally
		if ( $isolatetype eq "uks" ) {
		    $isolatetype = "se";
		}
	    }

	    $isolatetype_dominants{$isolatetype}++;
	}
    }
}

sub dhyper {
    my $M = shift; # white balls
    my $N = shift; # black balls
    my $m = shift; # white draws "success"
    my $n = shift; # black draws "fail"
	
    my $p = $M / ($M + $N);
    my $q = $N / ($M + $N);

    my $tries = $m + $n;
    
    if( $tries == 0 ) {
	$WARNING && output "\tProb undefined since count is 0.\n";
	next;
    }
	
    my $hypgeo;
    if( $m != 0 ) {
	$hypgeo = over($tries,$m) * ($p**$m) * ($q**$n) ;	
    } else {
	$WARNING && output "\tSucces count is 0; using fail count for hypergeometric prob.\n";
	$hypgeo = over($tries,$n) * ($p**$m) * ($q**$n) ;
    }

    my $m_over = "-";
    if ($m/$tries >= $p) {
	$m_over = "+";
    }
    
    if ($hypgeo eq "nan" or $hypgeo eq "inf") {
	$WARNING && output "\tOOPS - possibly too large numbers for dhyper: $hypgeo. Returning 1."; 
	$hypgeo = 1;
    }

    return ($hypgeo, $m_over);
}

sub over {
    my $a = shift;
    my $b = shift;
    
    if ($b == 0) { 
	return(undef);
    }

    my $facb = fac($b);
    defined($facb) || return(undef);

    my $ret = pfac($a-$b+1,$a) / $facb ;
}

sub pfac {
    my $y = shift;
    my $x = shift;

    if($y > $x or $x < 1) {
	return(undef);
    }
    

    my $ret = 1;
    foreach my $i ($y..$x ) {
	$ret *= $i;
    }

    return($ret);
}

sub fac {
    my $x = shift;
    
    my $ret = 1;
    if($x > 1) {
	foreach my $i ( 2..$x ) {
	    $ret *= $i;
	}
    }
	    
    if($x == 0) {
	$ret = undef;
    }
    
    return($ret);
}

# suggest segments

# first mark columns where id drops below treshold 

sub output {

    my $outstr = join("",@_);
    
    if ($LOG) {
	push @section_log, $outstr;
    } else {
	print $outstr;
    }
}

sub usage {
    print "USAGE: section_alignment.v2.pl\t<--in clustalw.aln>\n";
    print "\t\t\t\t[--dominancetab|-t dominance_table]\n";
    print "\t\t\t\t[--countmodel oneperisolate]\t\tcount one dominant per isolate\n";
    print "\t\t\t\t[--ultraconserved f ]\t0<f<1\t$ultraconserved_threshold\n";
    print "\t\t\t\t[--conserved f ]\t0<f<1\t$conserved_threshold\n";
    print "\t\t\t\t[--strictsection|-s]\tStrict section\t$strict_section\n";

    exit(1);
}
