#!/usr/bin/perl -w

use Bio::Matrix::PSM::IO;

# dump in rescore-using-R-readable format?

my $file = $ARGV[0];

my $dominance_tab = $ARGV[1];

my $sequence_file = $ARGV[2]; #fasta

#read_dominance

my %fractcount;

my $co_dominants = 0;
my $se_dominants = 0;

# to simplify reading, also get the fasta file?

use Bio::SeqIO;

my $fasta_file = Bio::SeqIO->new(-file => $sequence_file, '-format' => 'Fasta');

my %sequence; 
while ( my $seq = $fasta_file->next_seq() ) {
    $sequence{$seq->id} = $seq;
#    print "Sequence ",$seq->id," first 10 residues ",$seq->subseq(1,10),"\n";
}


open DOM, "<$dominance_tab";
while(my $r = <DOM>) {
    chomp $r;
    
    $r =~ /^(\S+)\s+(\S+)/;    
    $fractcount{$1} = $2;
   
    my $fractcount = $fractcount{$1};

    while ( $fractcount =~ m/co/g )  {
	$co_dominants++;
    }
	
    while ( $fractcount =~ m/cp/g ) {
	$co_dominants++;
    }

    while ( $fractcount =~ m/uks/g ) {
	$se_dominants++;
    }
    
    while ( $fractcount =~ m/se/g ) {
	$se_dominants++;
    }

}
close DOM;

print "From fraction-count: CO dominants: $co_dominants SE dominants: $se_dominants.\n";

print "Found 1 motifs.\n";

my $psmIO = new Bio::Matrix::PSM::IO(-format=>'meme', -file=>$file);

my $current_motif_id = "";
my $motifnr = 0;
my $current_motif_instance_lines ="";
my $motif_members = 0;

while (my $psm=$psmIO->next_psm) {
    my %psm_header=$psm->header;
    my $score=$psm_header{e_val};
    my $sites=$psm_header{sites};
    my $width=$psm_header{width};
    my $ic=$psm_header{IC};

    my $IUPAC=$psm->IUPAC;
    my $instances=$psm->instances;


    foreach my $instance (@{$instances}) {
	my $id=$instance->primary_id;
	
	$id =~ /^(\S+)\@(\S+)/;
	my $motif_id = $1;
	my $seq_id = $2;
	
	my $fractcount = $fractcount{$seq_id};

	if($motif_id ne $current_motif_id) {
	    # on each motif change 
	    
	    if ($motifnr != 0) {
		# on each change, but not on the first, dummy change from "" to motif1

		my $se_fraction = ($se_dominants_in_motif / ( $se_dominants_in_motif+$co_dominants_in_motif));
		
		my $se_over = ($se_fraction > 0.5) ? "+" : "-";
		my $dummy_pval = $score;

#		sprint "SE $se_dominants_in_mot\tCO\t$co_dominants\tSEfraction\t$se_fraction p $dummy_pval \n";

		print "Seed $current_motif_id - $motif_members members\n";
		print $current_motif_instance_lines;
		printf "SE %3d CO %3d SEfraction %.3f p %.4f %s\n", $se_dominants_in_motif, $co_dominants_in_motif, $se_fraction, $dummy_pval, $se_over;

		# reset motif counters
		$current_motif_instance_lines = "";
		$motif_members = 0; # no worries, will be increased again after this block
		print "\n"; 

	    } else {
		# on first motif 
#		print "Seed $motif_id - 0 members\n";
	    }

	    $motifnr++;

	    $co_dominants_in_motif = 0;
	    $se_dominants_in_motif = 0;
	}

	while ( $fractcount =~ m/co/g )  {
	    $co_dominants_in_motif++;
	}
	
	while ( $fractcount =~ m/cp/g ) {
	    $co_dominants_in_motif++;
	}

	while ( $fractcount =~ m/uks/g ) {
	    $se_dominants_in_motif++;
	}

	while ( $fractcount =~ m/se/g ) {
	    $se_dominants_in_motif++;
	}

	my $start_pos = $instance->start;
	my $end_pos = $start_pos + $width -1;
	my $score = $instance->score;
#	my $seq = $instance->seq;
	    
	my $subsequence = $sequence{$seq_id}->subseq($start_pos, $end_pos);
	$current_motif_instance_lines .= "$seq_id\t$start_pos\t$end_pos\t$score\t$subsequence\t".$fractcount."\n";
	
	$current_motif_id = $motif_id;
	$motif_members++;

	#Do something with the id
    }
}

# last motif dominant line..
#print "SE\t$se_dominants\tCO\t$co_dominants\n";
my $se_fraction = ($se_dominants_in_motif / ( $se_dominants_in_motif+$co_dominants_in_motif));
		
my $se_over = ($se_fraction > 0.5) ? "+" : "-";
my $dummy_pval = 1.1;

#		sprint "SE $se_dominants_in_mot\tCO\t$co_dominants\tSEfraction\t$se_fraction p $dummy_pval \n";

print "Seed $current_motif_id - $motif_members members\n";
print $current_motif_instance_lines;
printf "SE %3d CO %3d SEfraction %.3f p %.4f %s\n", $se_dominants_in_motif, $co_dominants_in_motif, $se_fraction, $dummy_pval, $se_over;

print "\nDone.\n";




