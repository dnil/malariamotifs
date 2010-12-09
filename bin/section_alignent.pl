#!/usr/bin/perl -w

use Bio::AlignIO;
use Bio::SearchIO;
use Bio::SeqIO;

use Bio::Tools::Run::StandAloneBlast;

sub read_fract_count;
sub find_peak;
sub check_remaining_hits_in_region;
sub run_blastpgp;
sub markup_domain;
sub check_domain_status;
sub get_unchecked_parts;

sub find_dominance;
sub get_best_hit_from_blast_gff_results;
sub exclude_certain_hits_from_dominance;
sub usage;

my @clargs = @ARGV;
my $DEBUG = 0;
my $WARNING = 1;

my $LOG = 1;

my $outfile = "";
my $infile = "";

my $fract_count = "";

my $dominance_exclusion_contig_ids_file = "";
my $dominance_exclusion_gbk_ids_file = "";
my $fractioncount_file = "";
my $dominance_model = "sum"; # the sum of primer contributions determines dominance order by default

# apparent blastpgp limitation (Or is that simply due to the "no hits found" bioperl issue?)
my $min_blastp_len = 6;

# default values
#my $min_init_w = 4;
#my $maxw = 40;
#my $ignore_section_treshold = 0;

my $dominant_cutoff = 5;

my $ultraconserved_threshold = 0.9;
my $conserved_threshold = 0.45;

my $tolerated_motif_overlap = 3;
my $hsp_neg_cutoff = 2;

my $count_one_dominant_per_isolate = 0;

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
#    } elsif($arg eq "--min_init_w") {	$arg = shift @clargs;
#	
#	if($arg eq "") {
#	    $WARNING && print STDERR "WARNING: no argument given to --min_init_w flag.\n";
#	} else {
#	    $min_init_w = $arg;
#	}
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
	    print "WARNING: unknown argument given to --countmodel.\n";	    
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
    } elsif($arg eq "--blastpgp_negatives") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no rows to delete given to --blastpgp_negatives flag.\n";
	} else {
	    $hsp_neg_cutoff = $arg;
	}
#    } elsif($arg eq "--maxw") {
#	$arg = shift @clargs;
#	
#	if($arg eq "") {
#	    $WARNING && print STDERR "WARNING: no rows to delete given to --maxw flag.\n";
#	} else {
#	    $maxw = $arg;
#	}
    } elsif(($arg eq "--fractioncount") || ($arg eq '-f')) {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no fraction count file given after --fractioncount flat.\n";
	} else {
	    $fractioncount_file = $arg;
	}

    } elsif($arg eq "--dominant_cutoff") {
	$dominant_cutoff = shift @clargs;
    } elsif($arg eq "--blasthitgff") {
	$blasthitgff_file = shift @clargs;
    } elsif($arg eq "--exclude_from_dominance") { # required for printing!
	$dominance_exclusion_gbk_ids_file = shift @clargs;
    } elsif($arg eq "--exclude_do_contigs") {
	$dominance_exclusion_contig_ids_file = shift @clargs;
    } elsif($arg eq "--dominance") {
	$dominance_model = shift @clargs;
	if( $dominance_model eq "sum" or $dominance_model eq "sum_or_af_top") {
	    $DEBUG && print STDERR "DEBUG: Using dominance model $dominance.\n";	    
	} else {
	    print STDERR "WARNING: Unrecognised dominance model $dominance. Using sum model instead.\n";
	    $dominance_model = "sum";
	}
    } elsif($arg eq "--debug") {
	$DEBUG = 1;
	$WARNING && print STDERR "DEBUG mode enabled.\n";
    } elsif($arg eq "--silent") {
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

# LOG execution parameters
$LOG && print "$0 ", join(" ", @ARGV), "\n";
$LOG && print "Run at ", scalar(localtime)," on ",$ENV{HOSTNAME}," in ",$ENV{PWD}," by ",$ENV{USER},"\n";

$LOG && print "min_blastp_len=$min_blastp_len\n";

#$LOG && print "min_init_w=$min_init_w\n";
#$LOG && print "maxw=$maxw\n";
#$LOG && print "ignore_section_treshold=$ignore_section_treshold\n";
$LOG && print "dominant_cutoff=$dominant_cutoff\n";
$LOG && print "countmodel oneperisolate=$count_one_dominant_per_isolate\n";
$LOG && print "ultraconserved_threshold=$ultraconserved_threshold\n";
$LOG && print "conserved_threshold=$conserved_threshold\n";

$LOG && print "tolerated_motif_overlap=$tolerated_motif_overlap\n";
$LOG && print "hsp_neg_cutoff=$hsp_neg_cutoff\n";

$LOG && print "blastpgp_e=$blastpgp_e\n";
$LOG && print "blastpgp_h=$blastpgp_h\n";
$LOG && print "blastpgp_f=$blastpgp_f\n";
$LOG && print "blastpgp_F=$blastpgp_F\n";
$LOG && print "blastpgp_b=$blastpgp_b\n";
$LOG && print "blastpgp_v=$blastpgp_v\n";

my @exclude_list;
my $exclude = 0;

if( $dominance_exclusion_gbk_ids_file ne "" ) {
    if( $blasthitgff_file eq "") {
	$WARNING && print "No blasthits file given so exclusion of hit ids from dominance computation will not be made.\n";    
	$exclude = 0;
    } else {
	$exclude = 1;
	open(EXCLUDE, "<$dominance_exclusion_gbk_ids_file");
	while (<EXCLUDE>) {
	    chomp;
	    if( $_ ne "" ) {
		print STDERR "Exclude ",$_, "\n";
		push @exclude_list, $_;
	    }
	}
	close EXCLUDE;
    }
}

if($dominance_exclusion_contig_ids_file ne "") {
    $exclude = 1;
    open(EXCLUDE_CONTIG, "<$dominance_exclusion_contig_ids_file");
    while(<EXCLUDE_CONTIG>) {
	chomp;
	if( $_ ne "" ) {
	    push @exclude_contig, $_;
	}
    }
    close EXCLUDE_CONTIG;
}

$ain = Bio::AlignIO->new('-file'=>$infile, '-format'=>'clustalw');
# $aout = Bio::AlignIO->new('-file'=>">$outfile", '-format'=>'clustalw');

my $aln = $ain->next_aln();

my $aln_start = 1; 
my $aln_end = $aln->length;

my $aln_height = $aln->no_sequences;

# read in the fraction count file, if given 
my %isolate_fraction;
my %isolate_contig_primer_fraction;
($fractioncount_file ne "") && (read_fract_count); 

my %blasthit;
($exclude == 1) && (get_best_hit_from_blast_gff_results);

my %exclude;
($exclude == 1) && (exclude_certain_hits_from_dominance);

my %dominant;
my %fractcount;
my $co_dominants;
my $se_dominants;
($fractioncount_file ne "") && (find_dominance);

$LOG && print "\nFrom fraction-count: CO dominants: $co_dominants SE dominants: $se_dominants.\n\n";

# prepare for later psiblast by making a blast-db

my $db_file = $infile.".db.fasta";
my $db_temp_seqio = Bio::SeqIO->new('-file'=>">$db_file", '-format'=>'Fasta');
foreach my $seq ($aln->each_seq) {
    my $seq_str = $seq->seq();
    $seq_str =~ s/-//g;
    $outseq = Bio::Seq->new( '-id' => $seq->id, '-seq' => $seq_str );
    $db_temp_seqio->write_seq($outseq);
}

system("formatdb -i $db_file -p T");

# clear the seed and matrix out name list files
open CLEAR, ">".$db_file.".sn";
close CLEAR;

open CLEAR, ">".$db_file.".pn";
close CLEAR;

# debug log
my @section_log;

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
	    $LOG && printf "Column: %3d  Res: %s  Count: %3d  Dom: %s\n", $i, $res, $count{$res}, $dominance_by_res{$res};
	} else {
	    $LOG && printf "Column: %3d  Res: %s  Count: %3d\n", $i, $res, $count{$res};
	}
    }

    push @res, $max_res;

    if($gap_count > $max_count) {
	push @gap, 1;
    } else {
	push @gap, 0;
    }

    $id_percent = $max_count / ($total_count - $gap_count);
#    $id_height = $max_count;

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

    $LOG && printf "\tGap: %1d  Gaps: %3d  Max count: %3d\n",$gap[$i-1],$gap_count, $max_count;
    $LOG && printf "\tid: %.3f\n",$id_percent,$max_count;
    $LOG && print "\n";
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

$LOG && print "res\t",$res,"\n";
$LOG && print "peak\t",$peak,"\n";
$LOG && print "idp\t",$idp,"\n";
$LOG && print "lidp\t",$lidp,"\n";
$LOG && print "gap\t",$gap,"\n";

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

$LOG && print "\nFind ultraconserved regions.\n\n";

my %motif_member_id;
my %motif_member_start; 
my %motif_member_end;

#my %motif_member_gapped_substr

my @used = (0) x scalar(@peak);
my $n_peaks = find_peak(@peak); # sets used..

my $found_motifs = scalar(keys %motif_member_id);
$LOG && print "Found $found_motifs different motifs in $n_peaks ultraconserved regions.\n";

my @non_gapped_idp;

for (my $i = 0; $i < scalar(@idp); $i++) {
    # drop the peaks that were already blasted from idp 
    # they are probably clean (residual mop-up is left for the last stage)
    $non_gapped_idp[$i] = ( ($idp[$i]> 0) && ($gap[$i]==0) && ($used[$i]==0) ? 1 : 0);
}

$LOG && print "\nFind ungapped high ID% regions, disregarding the previous peak motif parts.\n\n";

$LOG && print "res\t", $res,"\n";
$LOG && print "used\t",join("", @used),"\n";
$LOG && print "ng_idp\t", join("", @non_gapped_idp),"\n";

$n_peaks = find_peak(@non_gapped_idp);

my $current_motifs = scalar(keys %motif_member_id);
my $new_motifs = $current_motifs - $found_motifs;
$found_motifs = $current_motifs;

$LOG && print "Found $new_motifs different motifs in $n_peaks different ungapped high ID% regions.\n";

# now, directly look at the co-aligned, but unhit seqs?

my @gapped_idp;
for (my $i = 0; $i < scalar(@idp); $i++) {
    $gapped_idp[$i] = ((($idp[$i]>0) && ($gap[$i]>0)) ? 1 : 0);
}

$LOG && print "\nFind high ID% in gaps.\n\n";

$LOG && print "res\t", $res,"\n";
$LOG && print "g_idp\t", join("", @gapped_idp),"\n";
$LOG && print "gap\t", join("", @gap),"\n";

$n_peaks = find_peak(@gapped_idp);

$current_motifs = scalar(keys %motif_member_id);
$new_motifs = $current_motifs - $found_motifs;
$found_motifs = $current_motifs;

$LOG && print "Found $new_motifs different motifs in $n_peaks different gapped high ID% regions.\n";

#my @non_gapped_lidp;

#for (my $i = 0; $i < scalar(@lidp); $i++) {
#    $non_gapped_lidp[$i] = (($lidp[$i]> 0) && ($gap[$i]==0) && ($used[$i]==0));        
#}

# $DEBUG && print "\nFind ungapped high ID% regions, disregarding the previous peak motif parts.\n\n";

$LOG && print "\nLastly, find motifs in previously unhit parts.\n\n";

$DEBUG && print "res\t", $res,"\n";
$DEBUG && print "peak\t",$peak,"\n";
$DEBUG && print "ng_idp\t", join("", @non_gapped_idp),"\n";
$DEBUG && print "g_idp\t", join("", @gapped_idp),"\n";
$DEBUG && print "used\t", join("", @used),"\n";

my $used= join("",@used);
$used =~ s/1/9/g;
$used =~ s/0/1/g;
$used =~ s/9/0/g;

my @untouched = split(/ */,$used);
find_peak(@untouched);

$current_motifs = scalar(keys %motif_member_id);
$new_motifs = $current_motifs - $found_motifs;

$LOG && print "Found $new_motifs different motifs in $n_peaks different previously unhit regions.\n";

$DEBUG && print "\n\nStatus at end of loop:\n";
print "res\t", $res,"\n";

foreach my $id ( keys %status_ungapped ) {
    print $id, "\t", $status_ungapped{$id},"\n";
}

# construct new indicator where ungapped is 

# so, did we cover the peak and id intervening sections properly?
print @section_log;

print "\nMotifs\n";

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

    my $se_dominants_in_motif = 0;
    my $co_dominants_in_motif = 0;

    my $actual_members = scalar(@motif_id);
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
	    $WARNING && print STDERR "WARNING: no seq by the id ", $motif_id[$i],"\n";
	    $motif_bulkout[$motif_nr] .= $motif_id[$i]."\t".($motif_start[$i]+1)."\t".($motif_end[$i]+1);
	}

	if(defined($fractcount{$motif_id[$i]}) && $fractcount{$motif_id[$i]} ne "") {
	    
	    my $fractcount = $fractcount{$motif_id[$i]};

	    if($count_one_dominant_per_isolate == 1) {	
		while ( $fractcount =~ m/(co\d+)/g ) {
		    $current_count = $1;
		    if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {
			
		    } else {
			$counted_isolate{$current_count} = 1;
			$co_dominants_in_motif++;
		    }
		}
		
		while ( $fractcount =~ m/(cp\d+)/g ) {

		    $current_count = $1;
		    $current_count =~ s/cp/co/;
		    if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {
			
		    } else {
			$counted_isolate{$current_count} = 1;
			$co_dominants_in_motif++;
		    }
		}

		while ( $fractcount =~ m/(uks\d+)/g ) {
		    $current_count = $1;

		    if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {
			
		    } else {
			$counted_isolate{$current_count} = 1;
			$se_dominants_in_motif++;
		    }
		}
		
		while ( $fractcount =~ m/(se\d+)/g ) {
		    $current_count = $1;
		    if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {
			
		    } else {
			$counted_isolate{$current_count} = 1;
			$se_dominants_in_motif++;
		    }
		}

	    } else {

		while ( $fractcount =~ m/co\d+/g )  {
		    $co_dominants_in_motif++;
		}

		while ( $fractcount =~ m/cp\d+/g ) {
		    $co_dominants_in_motif++;
		}

		while ( $fractcount =~ m/uks\d+/g ) {
		    $se_dominants_in_motif++;
		}
		
		while ( $fractcount =~ m/se\d+/g ) {
		    $se_dominants_in_motif++;
		}
	    }

	    $motif_bulkout[$motif_nr] .= "\t".$fractcount;

	}
	$motif_bulkout[$motif_nr] .= "\n";
    }

    if($se_dominants_in_motif + $co_dominants_in_motif > 0 ) {
	$DEBUG && print "DEBUG: Calling dhyper for seed $seed_id ($actual_members actual members) with M=$se_dominants, N=$co_dominants, m=$se_dominants_in_motif, n=$co_dominants_in_motif.\n";
	my ($pval, $se_over) = dhyper($se_dominants, $co_dominants, $se_dominants_in_motif, $co_dominants_in_motif);
	$motif_pval[$motif_nr] = $pval;
#	$motif_seover[$motif_nr] = $se_over;
	$motif_bulkout[$motif_nr] .= sprintf "SE %3d CO %3d SEfraction %.3f p %.4f %s\n", $se_dominants_in_motif, $co_dominants_in_motif, ($se_dominants_in_motif / ( $se_dominants_in_motif+$co_dominants_in_motif)), $pval, $se_over;
#	exec("echo |R")
    }
    $motif_nr++;
}

# sort
print "Found $motif_nr motifs.\n";

my @order;
if ($fractioncount_file ne "") {
    @order = sort {$motif_pval[$a] <=> $motif_pval[$b] } (0..($motif_nr-1));   
} else {
    @order = 0..($motif_nr-1);
}

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
	$DEBUG && print "Sequence ",$seed_seq->id," has no perfect match (at start)..\n";
	return ("", 0,0,0,0);
    }
    my $peak_seq_start = $peak_seq_loc_start->start; # 1-based inclusive
    
    my $peak_seq_loc_end = $seed_seq->location_from_column($peak_end);
    if (!defined($peak_seq_loc_end) || $peak_seq_loc_end->location_type() eq 'IN-BETWEEN') {
	$DEBUG && print "Sequence ",$seed_seq->id," has no perfect match (at end)..\n";
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
	    $DEBUG && print "Examining pos $i; not at peak at the moment.\n";
	    if($indicator[$i-1] == 1) {
		$peak_start = $i;
		$peak_end = $i;
	    }
	} else {
	    if($indicator[$i-1] == 1) {
		$peak_end = $i;
	    } elsif(($indicator[$i-1] == 0) || ($i == $aln_end)) {
		# leaving peak

		if( $peak_start == $peak_end ) {
		    $DEBUG && print "Ignoring single peak residue at pos $peak_start.\n";

		    $peak_start= -1;
		    $peak_end= -1;
		    next;
		}

#		$LOG && print "Found indicator peak between $peak_start and $peak_end.\n";

		push @section_log, "Indicator peak: $peak_start - $peak_end.\n";
		$n_peaks++;
#		foreach my $setpos ( $peak_start..$peak_end ) {
#		    $touched[$setpos-1] = 1;
#		}
		
		# 1. pick a peak-representative read
		
		my $seed_seq;
		my $found_representative = 0;
		
		my $peak_string = substr($res, $peak_start-1, $peak_end-$peak_start+1);
		
		$DEBUG && print "Looking for representative for $peak_string.. ($peak_start - $peak_end)\n";
		for(my $j = 1; $j <= $aln_height; $j++) {

		    $seed_seq = $aln->get_seq_by_pos($j);
		    		    
		    my ($rep_cand_str, $rep_cand_len, $peak_seq_start, $peak_seq_end, $peak_len) = subseq_from_columns($seed_seq, $peak_start, $peak_end);
		    if ($rep_cand_str eq "") {
			next; # failed
		    }
		    $rep_cand_substr = substr($rep_cand_str, $peak_seq_start-1, $peak_len); # substr takes 0-based coord, peak coord is 1-based

		    # check status!
		    my $tried = check_domain_status($seed_seq->id, $peak_seq_start, $peak_seq_end); # call with 1-based incl.

		    if($tried > 0) { 
			$WARNING && print STDERR "Candidate contig ", $seed_seq->id, " $peak_seq_start - $peak_seq_end has status $tried. Not a rep.\n";
			next;
		    }

		    if ($rep_cand_len < $min_blastp_len) {
			$WARNING && print STDERR "Candidate contig ", $seed_seq->id, " is shorter $rep_cand_str than min blastp query ($min_blastp_len).\n";
			next;
		    }

		    $DEBUG && print "Sequence ",$seed_seq->id," $peak_seq_start - $peak_seq_end : $rep_cand_substr";
		    
		    if( $rep_cand_substr eq $peak_string) {
			$DEBUG && print "..match!\n";

			$found_representative = 1;

			# Extend hit if too short for blast
			if($peak_len < $min_blastp_len) {

			    ($peak_len,$peak_seq_start, $peak_seq_end) = extend_hit_region($peak_len, $peak_seq_start, $peak_seq_end, $rep_cand_str, $rep_cand_len);
	    
			    $rep_cand_substr = substr($rep_cand_str, $peak_seq_start-1, $peak_seq_end-$peak_seq_start+1);
			    $peak_len = $peak_seq_end-$peak_seq_start+1;

			    $WARNING && print "Peak ($peak_string) too short for blastp; using ", $seed_seq->id,
			    " $peak_seq_start - $peak_seq_end ($rep_cand_substr) as a representative.\n";
			}

			# mark query as tested -- really not needed.. Should hopefully be hit by its own representative.. =)
			my $save = markup_domain($seed_seq->id, $peak_seq_start, $peak_seq_end, 1);

			# also mark this region in consensus as "used"; so far only used to remove "peaks" from "idp"-regions
			@used[($peak_start-1)..($peak_end-1)] = (1) x ($peak_end - $peak_start+1);
			
			# construct representative for blasting
			my $old_id = $seed_seq->id;
			$seed_seq = Bio::Seq->new( '-id' => $old_id."_".$peak_seq_start."_".$peak_seq_end, '-seq' => $rep_cand_substr);
			if($save == 1) {
			    save_motif($seed_seq->id, $old_id, $peak_seq_start, $peak_seq_end);
			    $DEBUG && print "Saving motif ",$seed_seq->id,", $old_id, $peak_seq_start, $peak_seq_end [find_peak]\n";
			}

			last;
		    } else {
			$DEBUG && print "..mismatch!\n";
		    }
		}
		
		if ($found_representative == 0) {
#		    $WARNING && print "Could not find exact representative for peak at columns $peak_start - $peak_end ($peak_string). Proceeding to next motif instead.\n";
		    $WARNING && print "Could not find exact representative for peak at columns $peak_start - $peak_end ($peak_string). Proceeding to by row motif search instead.\n";
		    # next;		    
		    
		} else {		    
		    # 2. run blastpgp
		    		   
		    run_blastpgp($seed_seq, $db_file);
		}

		$DEBUG && print "\nCheck remaining hits in columns $peak_start - $peak_end..\n";
		check_remaining_hits_in_region($peak_start, $peak_end);
		$DEBUG && print "\nLeft checking remaining hits $peak_start - $peak_end..\n";

		$peak_start= -1;
		$peak_end= -1;

		# then, immediately check any remaining non-hits as queries, until all have been tried
		
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
	$peak_seq_end = $rep_cand_len;
    } else {
	$peak_seq_end += $needed_extension_last_half;
    }
    
    # update candidate substrin
    $rep_cand_substr = substr($rep_cand_str, $peak_seq_start-1, $peak_seq_end-$peak_seq_start+1);
    $peak_len = $peak_seq_end-$peak_seq_start+1;

    return ($peak_len,$peak_seq_start, $peak_seq_end);
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
      my $seed_seq_reg_str = substr($seed_seq_str, $seed_seq_start-1, $seed_seq_end-$seed_seq_start+1);
      # check status!
      my $tried = check_domain_status($seed_seq->id, $seed_seq_start, $seed_seq_end);
	
      if($tried > 0) {
	  $WARNING && print "Candidate contig ", $seed_seq->id, " $seed_seq_start - $seed_seq_end has status $tried. No need for additional search!\n";
	  next ROW; # get next row instead
      }

      $DEBUG && print "\nCheck remaining hits says sequence ",$seed_seq->id," $seed_seq_start - $seed_seq_end : $seed_seq_reg_str was untouched.\n";

      # Extend hit if too short for blast
	
      if ($seed_seq_len < $min_blastp_len) {
	  $WARNING && print "Candidate ", $seed_seq->id, " is shorter $seed_seq_str than min blastp query ($min_blastp_len).\n";
	  next;
      }

	if($peak_len < $min_blastp_len) {
	    $DEBUG && print "Attempting to extend...\n";
	    ($peak_len,$seed_seq_start, $seed_seq_end) = extend_hit_region($peak_len, $seed_seq_start, $seed_seq_end, $seed_seq_str, $seed_seq_len);
	    $seed_seq_reg_str = substr($seed_seq_str, $seed_seq_start-1, $seed_seq_end-$seed_seq_start+1);

	    $WARNING && print "Peak too short (now $peak_len) for blastp; using ", $seed_seq->id,
	    " $seed_seq_start - $seed_seq_end ($seed_seq_reg_str) as a representative.\n";
	}

	# mark query as tested 
	my $save = markup_domain($seed_seq->id, $seed_seq_start, $seed_seq_end, 1);

	# no need to save status -- really not needed.. Should hopefully be hit by its own representative.. =)

	my $old_id = $seed_seq->id;
	$seed_seq = Bio::Seq->new( '-id' => $old_id."_".$seed_seq_start."_".$seed_seq_end, '-seq' => $seed_seq_reg_str);

	if($save == 1) {
	    save_motif($seed_seq->id, $old_id, $seed_seq_start, $seed_seq_end);
	    $DEBUG && print "Saving motif ",$seed_seq->id,", $old_id, $seed_seq_start, $seed_seq_end [from check remaining hits]\n";
	}

	run_blastpgp($seed_seq, $db_file);

    }
}

sub run_blastpgp {

    my $seq = shift;
    my $database = shift; # first, make sure there is a database..
    my $seq_name = $seq->id;
    my $outfile_basename = $database."_searched_by_".$seq_name;
    
    # write the sequence file... filename .ceq
    my $ceq_file = $outfile_basename.".ceq";

    open CEQ, ">$ceq_file";
    print CEQ ">$seq_name\n";
    print CEQ $seq->seq(),"\n";
    close CEQ;

    # add ceq-file name to list (kind of slow to use append; the alternative is to keep an open filehandle from caller)
    open SNLIST, ">>$database.sn";
    print SNLIST $ceq_file."\n";
    close SNLIST;

    my $chk_file = $outfile_basename.".blpgp.chk"; 

    # add the to be produced chk-file name to list (kind of slow to use append; the alternative is to keep an open filehandle from caller and pass it, or use a global glob)
    open PNLIST, ">>$database.pn";
    print PNLIST $chk_file, "\n";
    close PNLIST;

    @params = ('database' => $database,'outfile' => $outfile_basename.".blpgp",
	       'e' => $blastpgp_e, 'h'=>$blastpgp_h, 'f'=>$blastpgp_f, 'F'=>$blastpgp_F, 'b'=>$blastpgp_b, 'v'=>$blastpgp_v,
	       # testade e=1000 med dåliga resultat..
#	       'm'=>6,
	       'C'=>$chk_file,
	       '_READMETHOD'=> 'BPlite'
	       );
# m=6 is a flag to get aln working according to BPlite perldoc
    $DEBUG && print "Blastpgp of $seq_name (".($seq->seq).") vs $database..\n";

    my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
    $factory->j(10);
    my $report = $factory->blastpgp($seq);
    
    if( $report->{'REPORT_DONE'} == 1 ) {
	$WARNING && print "OOPS: Report reported done before calling number_of_iterations.\n";
	$WARNING && print "This probably means there were no hits below e-val treshold. Dropping current query.\n";
	return;
    }

    my $total_iterations = $report->number_of_iterations;
    my $last_iteration = $report->round($total_iterations);

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

    my $hsp_pos_cutoff = $last_iteration->qlength - $hsp_neg_cutoff;
    
    my @markupqueue;
    while (my $sbjct = $last_iteration->nextSbjct) {

	$DEBUG && print "Iterated blast hit ", $sbjct->name,".\n";
	
	while (my $hsp = $sbjct->nextHSP) {
	    if($hsp->positive() >= $hsp_pos_cutoff) {
		$DEBUG && printf "Motif member %10s  %3d - %3d : %s\n", $hsp->hit->seq_id, $hsp->hit->start, $hsp->hit->end, $hsp->sbjctSeq;		

		push @markupqueue, [$hsp->hit->seq_id, $hsp->hit->start, $hsp->hit->end, 2];
		
	    } else {
		$DEBUG && print "Dropping HSP since positives ",$hsp->positive," less than $hsp_pos_cutoff (",$hsp->sbjctSeq,") .\n";
	    }
	}
    }
    if (@markupqueue > 1) {
	# mark subject hit!
	while (my $markme = shift @markupqueue) {
	    my $save = markup_domain(@$markme);
	    if($save == 1) {
		save_motif($seq_name, @$markme);
		$DEBUG && print "Saving motif $seq_name ",join(" ",@$markme)," from run_blastpgp for seed ".$$markme[0]."\n";
	    }
	}
    } else {
	# drop this alignment/motif
	$WARNING && print "Number of subjects hit < 2, so dropping current query.\n";
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
	    $WARNING && print "WARNING: Extensive double handling of $id $start - $end: was $prev is ",($type x $hit_len),". Ignoring match.\n";
	    return 0;
	}
    }
    # actuall replace

    $DEBUG && print "setting status for $id $start - $end to ",($type x $hit_len),".\n";
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
	    $WARNING && print "WARNING: Previously touched $id $start - $end: $status. Ignoring potential seed.\n";
	    return(1);
	}
	$WARNING && print "WARNING: Previously touched $id $start - $end: $status. Using anyway.\n";
	return(0);
    } else {
	$DEBUG && print "DEBUG Checking status for $id $start - $end: $status. Untouched.\n";
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

sub read_fract_count {

    $DEBUG && print "Attempting to open fraction count file $fractioncount_file.\n";
    open FRACTCOUNT,"<$fractioncount_file ";
    while(<FRACTCOUNT>) {
	chomp;
	my @row = split /\s+/, $_;

	my $isolate = $row[2];
	my $contig = $row[0];
#    my $isolate_fraction = $row[5];
	my $isolate_fraction = $row[7];
	${$isolate_fraction{$isolate}}->{$contig} = $isolate_fraction;

	if ( $dominance_model ne "sum" ) {
	    my $primer = $row[5];
	    my $isolate_primer_fraction = $row[8];
    
	    ${${$isolate_contig_primer_fraction{$isolate}}->{$contig}}->{$primer} = $isolate_primer_fraction;
	}

#    $DEBUG && print "DEBUG: $isolate $contig $isolate_fraction\n";
    }   
}


sub find_dominance {

    # reset counters (global scope)
    $co_dominants = 0;
    $se_dominants = 0;

    foreach my $isolate (keys %isolate_fraction) {
	
	my (@contigs) = keys %${$isolate_fraction{$isolate}};
	my $maxcount = ( @contigs < $dominant_cutoff ) ? @contigs : $dominant_cutoff;

	if ( $exclude == 1 ) {
	    my @include_contigs;
	    foreach my $contig (@contigs) {
		if ( defined( $exclude{$contig} ) ) {
		    if( $exclude{$contig} == 1 ) {
			$DEBUG && print STDERR "Excluded contig $contig.\n";
		    } elsif ($exclude{$contig} == 0) {
			push @include_contigs, $contig;
		    }
		} else {
		    # if contig didn't get a blast hit, or elsewise was marked excluded, mark it non-excluded =)
		    $exclude{$contig} = 0; 
		    # AND includet, you nipwit
		    push @include_contigs, $contig;
		}
	    }
	    @contigs = @include_contigs;
	}
	
	if( $dominance_model eq "sum" ) { 
	    
	    @contigs = sort { ${$isolate_fraction{$isolate}}->{$b} <=> ${$isolate_fraction{$isolate}}->{$a}} @contigs;
	    
	    $DEBUG && print "$maxcount dominants to process for $isolate (out of ",scalar(@contigs)," contigs).\n";

#    $DEBUG && print "DEBUG: $isolate maxcount $maxcount\n";
	
	    my %counted_isolate;   # once for the entire contig list - count nr of isolates

	    for( my $i = 0; $i < $maxcount ; $i++ ) {
		$dominant{$contigs[$i]} = 1;
		if( exists($fractcount{$contigs[$i]}) ) {
		    # this is a (sub)dominant fraction
		    $fractcount{$contigs[$i]} .= "_".$isolate.".".($i+1)."_".${$isolate_fraction{$isolate}}->{$contigs[$i]};
#	    $DEBUG && print "DEBUG: ".$contigs[$i]." ".$fractcount{$contigs[$i]}."\n";
		} else {
		    $fractcount{$contigs[$i]} = $isolate.".".($i+1)."_".${$isolate_fraction{$isolate}}->{$contigs[$i]};		
#	    $DEBUG && print "DEBUG: ".$contigs[$i]." ".$fractcount{$contigs[$i]}."\n";
		}

	    }      
	} elsif ( $dominance_model eq "sum_or_af_top" ) {
	
	    my $afdominant_ok = 0;
	    my $afdominant = "";  

	    @contigs = sort { ${$isolate_fraction{$isolate}}->{$b} <=> ${$isolate_fraction{$isolate}}->{$a}} @contigs;
	
	    my %contigs_primer;

	    foreach my $contig ( @contigs ) {
		foreach my $primer ( "af","nfbr","ndbl" ) { # ehrr, why not obtain keys from the hash?
		    if ( defined(${${$isolate_contig_primer_fraction{$isolate}}->{$contig}}->{$primer}) ) {
			push @{$contigs_primer{$primer}}, $contig;
		    }
		}
	    }

	    if( defined($contigs_primer{af}) && @{$contigs_primer{af}} > 0 ) {
	    
		@{$contigs_primer{af}} = sort { ${${$isolate_contig_primer_fraction{$isolate}}->{$b}}->{af} <=> ${${$isolate_contig_primer_fraction{$isolate}}->{$a}}->{af} } @{$contigs_primer{af}};	
		
		$afdominant = ${$contigs_primer{af}}[0];
	
		$afdominant_ok = 0;
		for (my $i = 0 ; $i < $maxcount; $i++ ) { 
		    if ( $contigs[$i] eq $afdominant ) { 
			# afdominant is among the top N - no worries
			$afdominant_ok = 1;
			$DEBUG && print "DEBUG: Isolate $isolate has the top af dominant among the top $maxcount sum-dominants.\n";
		    }
		}
	    
	    } else {
		$afdominant_ok = 1;
		$DEBUG && print "DEBUG: Isolate $isolate has no af primer contigs at all. Using $maxcount sum-dominants from ndbl/nfbr only.\n";
	    }
	    
	    if ( $afdominant_ok == 1) {

		for (my $i = 0 ; $i < $maxcount; $i++ ) { 
		    if( exists($fractcount{$contigs[$i]}) ) {
			$fractcount{$contigs[$i]} .= "_".$isolate.".".($i+1)."_".${$isolate_fraction{$isolate}}->{$contigs[$i]};

			$DEBUG && print "DEBUG: ".$contigs[$i]." ".$fractcount{$contigs[$i]};
			if ($contigs[$i] eq $afdominant) { 
			    $DEBUG && print " *afdominant*\n";
			    
			} else {
			    $DEBUG && print "\n";
			}
		    } else {
			$fractcount{$contigs[$i]} = $isolate.".".($i+1)."_".${$isolate_fraction{$isolate}}->{$contigs[$i]};
			$DEBUG && print "DEBUG: ".$contigs[$i]." ".$fractcount{$contigs[$i]};
			if ($contigs[$i] eq $afdominant) {
			    $DEBUG && print " *afdominant*\n";
			} else {
			    $DEBUG && print "\n";
			}
		    }
		}
	    } else { 
		if( exists( $fractcount{$afdominant}) ) {
		    $fractcount{$afdominant} .= "_".$isolate.".1"."_".${$isolate_fraction{$isolate}}->{$afdominant} ;
		} else {
		    $fractcount{$afdominant} = $isolate.".1"."_".${$isolate_fraction{$isolate}}->{$afdominant} ;
		}
		for (my $i = 0 ; $i < $maxcount -1; $i++ ) { 
		    if( exists($fractcount{$contigs[$i]}) ) { 
			$fractcount{$contigs[$i]} .= "_".$isolate.".".($i+2)."_".${$isolate_fraction{$isolate}}->{$contigs[$i]} ;
		    } else {
			$fractcount{$contigs[$i]} = $isolate.".".($i+2)."_".${$isolate_fraction{$isolate}}->{$contigs[$i]} ;
		    }
		}
	    }
	}



	# NOTE: change for one-per-isolate count tactic..
	
	# and, for the final P-val, count nr of co (count each isololate only once) over nr of se isolates +  nr of co isolates ?
	# problem is, with only 6 dominants per isolate, we could run out of them (equal balls in urn)

            # but probably ok: the logic being that when we pull, say, two contigs dominant in co01, we count those as one (they should perhaps
	    # have been joined in assembly)



# The count of one-per-isolate balls is a bit more dodgy. I currently think the most correct would be to first assigning all the motifs, then count them through, with 
# any same-isolate-dominant-in-several-contigs-in-motif counting only once both in the totals and in the per motif tally. 
# No. Trouble is, that gives a count for each contig each time it appears in any motif, which may be several. So the initial bag of dominants will appear much larger than it is in reality.
# It's best to concede that the hypergeo-model is not applicable to the problem. 

# 	    if($count_one_dominant_per_isolate == 1) {	
# 		my $current_count;
# 		while ( $fractcount =~ m/(co\d+)/g ) {
# 		    $current_count = $1;
# 		    if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {
			
# 		    } else {
# 			$counted_isolate{$current_count} = 1;
# 			$co_dominants++;
# 		    }
# 		}
		
# 		while ( $fractcount =~ m/(cp\d+)/g ) {
# 		    $current_count = $1;
# 		    $current_count =~ s/cp/co/;
# 		    if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {
			
# 		    } else {
# 			$counted_isolate{$current_count} = 1;
# 			$co_dominants++;
# 		    }
# 		}
		
# 		while ( $fractcount =~ m/(uks\d+)/g ) {
# 		    $current_count = $1;
# 		    if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {
			
# 		    } else {
# 			$counted_isolate{$current_count} = 1;
# 			$se_dominants++;
# 		    }
# 		}
		
# 		while ( $fractcount =~ m/(se\d+)/g ) {
# 		    $current_count = $1;
# 		    if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {
			
# 		    } else {
# 			$counted_isolate{$current_count} = 1;
# 			$se_dominants++;
# 		    }
# 		}
		
# 	    } else {

    }

# 101011 Updating for Cameroon sequences
# "CO" cbm, "SE" cbsa, cbcs
    
    foreach my $contig (keys %fractcount) {
	my $fractcount = $fractcount{$contig};
	while ( $fractcount =~ m/cbm/g )  {
	    $co_dominants++;
	}
	
	while ( $fractcount =~ m/cbsa/g ) {
	    $se_dominants++;
	}

	while ( $fractcount =~ m/cbcs/g ) {
	    $se_dominants++;
	}
    }
   
}



sub get_best_hit_from_blast_gff_results {

    open BLASTGFF, "$blasthitgff_file";

    my %lowestexpect;
    
    while( my $row = <BLASTGFF> ) {

	my @col = split /\t+/, $row; 

	my $contigname = $col[0];
	my $querystart = $col[3]; 
	my $queryend = $col[4];
	my $colelements = scalar(@col)-1;
	my $comment = join "\t", @col[8..$colelements];

	my @comment = split /;+/, $comment;
	
	#   print STDERR $comment[0]."\n";
	
	# HIGHLY MALARIA DATASET-SPECIFIC!
	my $subjectname;
	if ($comment[0] =~ m/b\|([\w\d\.]+)\|/) {
	    $subjectname = $1;
	} else {
	    ($subjectname) = ($comment[0] =~ m/([\w\d\_\.]+)\|mRNA\|/);
	}

	my ($expect) = ($comment[1] =~ m/Expect = ([\deE\.\-]+)/);
	my ($querycoverage) = ($comment[1] =~ m/Querycoverage = ([\d\.]+)/);
	my ($idfreq) = ($comment[1] =~ m/Ids = \d+\/\d+ \(([\d\.]+)\)/);

	if( !defined($lowestexpect{$contigname}) || ($expect < $lowestexpect{$contigname}) ) { 
	    $blasthit{$contigname} = "_".$subjectname."_$querystart+".$queryend."_$querycoverage"."_".$idfreq;
	    $lowestexpect{$contigname} = $expect;
	    $DEBUG && print STDERR $contigname." ".$blasthit{$contigname}." ".$lowestexpect{$contigname}."\n";
	}

#    $DEBUG && print STDERR $contigname." ".$blasthit{$contigname}." ".$lowestexpect{$contigname}."\n";    
    }
}

sub exclude_certain_hits_from_dominance {

    $DEBUG && print STDERR "checking excludes for ", scalar(keys %blasthit), " hit contigs.\n";
    
    # check exclusion id list fo to 
    foreach my $contig (keys %blasthit) {

	$exclude{$contig} = 0;

	foreach my $exclude_id (@exclude_list) {
	    if ($blasthit{$contig} =~ m/\_$exclude_id\_/) {
		$DEBUG && print STDERR "Excluding hit for $contig to excluded id $exclude_id.\n";
		$exclude{$contig} = 1;
	    }
	}
    }

    if( $dominance_exclusion_gbk_ids_file ne "" ) {	    
	foreach my $exclude_contig (@exclude_contig ) {
	    $exclude{$exclude_contig} = 1;
	    $DEBUG && print STDERR "Excluding $exclude_contig based on contig exclude list entry.\n";
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
	$WARNING && print "\tProb undefined since count is 0.\n";
	next;
    }
	
    my $hypgeo;
    if( $m != 0 ) {
	$hypgeo = over($tries,$m) *  $p**$m * $q**$n ;
    } else {
	$WARNING && print "\tSucces count is 0; using fail count for hypergeometric prob.\n";
	$hypgeo = over($tries,$n) *  $p**$m * $q**$n ;
    }

    my $m_over = "-";
    if ($m/$tries >= $p) {
	$m_over = "+";
    }
    
    if ($hypgeo eq "nan" or $hypgeo eq "inf") {
	$WARNING && print "\tOOPS - possibly too large numbers for dhyper: $hypgeo. Returning 1."; 
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

sub usage {
    print "USAGE: section_alignment.pl\t<--in clustalw.aln>\n";
    print "\t\t\t\t[--fractioncount|-f frac_sort]\n";
    print "\t\t\t\t[--dominant_cutoff int]\t\t$dominant_cutoff\n";
    print "\t\t\t\t[--countmodel oneperisolate]\t\tcount one dominant per isolate\n";
    print "\t\t\t\t[--blasthitgff my.bln.gff]\n";
    print "\t\t\t\t[--exclude_from_dominance exclude_gbk_id_file]\n";
    print "\t\t\t\t[--ultraconserved f ]\t0<f<1\t$ultraconserved_threshold\n";
    print "\t\t\t\t[--conserved f ]\t0<f<1\t$conserved_threshold\n";

    exit(1);
}
