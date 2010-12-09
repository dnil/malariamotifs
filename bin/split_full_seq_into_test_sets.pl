#!/usr/bin/perl -w

my $DEBUG = 0;
my $WARNING = 1;

my $fastafile = "";
my $set1_fract = 2/3; 
my $qualfile = "";
use POSIX;

while (my $arg = shift @ARGV) {

    if ($arg =~ /^-/) {	
	if($arg eq '-D' or $arg eq '--debug') { 
	    $DEBUG = 1;
	} elsif ($arg eq '-n') {
	    $arg = shift @ARGV;
	    $N = $arg;
	    $set1_fract = 1/$N; # test
	} elsif ($arg eq '-q') {
	    $arg = shift @ARGV;
	    $qualfile = $arg;
	}
    } else {
	$fastafile = $arg;
    }
}

if($fastafile eq "") {
    print STDERR "No fasta file to split was given. Nothing to do!\n";
    exit 1;
}

if ($qualfile eq "") {
    print STDERR "FATAL: no qual file name given! If you're going to phrap this, you'll need to make sure that the out seq and qual files have entries in the exact same order, so you'd better specify one by giving -q qualfilename. Bailing out.\n";
    exit 1;
}

open FASTAFILE, $fastafile or die "Could not open $fastafile.\n";
open QUALFILE, $qualfile; 

my $set2_fract = 1 - $set1_fract; #train

my $printentry = 0;

my %seq; 
my %header;

my %patients;
my %readtype;  # Counter for the number of READS of each isolate type (#se/#uks/#mi)

#my %patienttype;  # Counter for the number of PATIENTS of each isolate type

my $current_read = "";
my $patient;
my ($type, $nr);

my %typesort; # hash of arrays of patients, keyed on type
my %patientsort; # hash of arrays of reads, keyed on patient

while( my $row = <FASTAFILE> ) {
    chomp $row;
    if ( $row =~ m/^>/ ) {
	$DEBUG && print STDERR "Checking $row.\n";
	if(($current_read) = ($row =~m/^>(\S+)/) ) {
	   $DEBUG && print $current_read,"\n";
	   
	   ($type, $nr) = ($row =~ /[rf]([a-zA-Z]{2,3})(\d{1,2})\w+\d+\.\w+\s+/);
	   $type = lc $type;
	   $nr = lc $nr;
	   $patient = $type.$nr;
	   
	   if( !defined( $readtype{$type} ) ) {
	       $readtype{$type} = 1;
	   } else {
	       $readtype{$type}++; # number of READs of type type
	   }

#	   if( !defined( $patienttype{$type} ) ) {
#	       $patienttype{$type} = 1; 
#	   } else {
#	       $patienttype{$type}++; # number of READs of type type
#	   }
	   
	   if ( !defined( @{$patientsort{$patient}} ) ) {
	       # first time encountering a read from this patient
	       # add patient to typesort
	       push @{$typesort{$type}}, $patient;	   
	       @{$patientsort{$patient}} = ();
	   }

	   push @{$patientsort{$patient}}, $current_read;

	   $header{$current_read} = $row;
	   $seq{$current_read} = "";

	} else {
	    $WARNING && print "WARNING: Fasta header row did not contain the magic read name format - ignoring $row.\n";
	}

    } else {
	if ($current_read ne "") {
	    $seq{$current_read} .= $row;
	    $DEBUG && print "added ", length($row), "bases to sequence $current_read.\n";
	}
    }
}


while( my $row = <QUALFILE> ) {
    chomp $row;
    if ( $row =~ m/^>/ ) {
	$DEBUG && print STDERR "Checking $row.\n";
	if(($current_read) = ($row =~m/^>(\S+)/) ) {
	   $DEBUG && print $current_read,"\n";
	   
#	   ($type, $nr) = ($row =~ /[rf]([a-zA-Z]{2,3})(\d{1,2})\w+\d+\.\w+\s+/);
#	   $type = lc $type;
#	   $nr = lc $nr;
#	   $patient = $type.$nr;
	   
	   $qual{$current_read} = "";

	} else {
	    $WARNING && print "WARNING: QUAL file fasta header row did not contain the magic name - ignoring $row.\n";
	}
    } else {
	if ($current_read ne "") {
	    $qual{$current_read} .= $row."\n";
	    $DEBUG && print "added ", length($row), "chars to qual string for $current_read.\n";
	}
    }
}
# now, first check how many patients we have in each category co/se/uks - ok done while reading file.

# then count how many of each an Nth is. 

# also, for randomisation we want one list for each type, of size $readtype{$type}
# and we don't want independent pickings really; we want to have each patient only once in the test set.. 

# 1. randomize patient order

my %patient_list; 

foreach $type (keys %readtype) {

    # how large fractions?
    my $nr_of_patients_of_this_type = @{$typesort{$type}};
    $set1_patients{$type} = floor( $set1_fract * $nr_of_patients_of_this_type ); # test
    $set2_patients{$type} = $nr_of_patients_of_this_type - $set1_patients{$type}; #train

#    $patient_counter{$type} = 0;

    my @patient_ids_for_type = @{$typesort{$type}};

    for (my $i = 0; $i < $nr_of_patients_of_this_type; $i++) {
	my $pick_next = int(rand($nr_of_patients_of_this_type - $i)); 

	push @{$patient_list{$type}}, $patient_ids_for_type[$pick_next];
	splice @patient_ids_for_type, $pick_next, 1;
    }
}

# 2a. check sets-offset for any part of set1 to pick first? 
# 2b. pick a set2 (the small fraction!). increment some offset so that t

open LOG, ">$fastafile.setsplit.log";

for (my $sets=0; $sets < $N ; $sets++) {

    open PFH, ">$fastafile.test.$sets.set";
    open QFH, ">$fastafile.train.$sets.set";

    open PQFH, ">$fastafile.test.$sets.set.qual";
    open QQFH, ">$fastafile.train.$sets.set.qual";

    print LOG "\n\nset $sets\n\n";

    foreach $type (keys %readtype) {

	print LOG "isolate type $type should have ".$set1_patients{$type}." test and ".$set2_patients{$type}." training set members.\n";
       
	for (my $patient_counter = 0;  $patient_counter < @{$patient_list{$type}} ; $patient_counter++ ) {

	    my $patient = ${$patient_list{$type}}[$patient_counter];
	    
	    if($patient_counter >= $sets * $set1_patients{$type} && $patient_counter < ($sets+1) * $set1_patients{$type}) {
		# This is a member of the current set1, test set, p-set!
		print LOG "entry set $sets, test-set member $patient: ".scalar(@{$patientsort{$patient}})." reads.\n";
		$outfh = "PFH"; 
		$qualoutfh = "PQFH";
	    } else {
		print LOG "set $sets, training-set member $patient: ".scalar(@{$patientsort{$patient}})." reads.\n";
		$outfh = "QFH";
		$qualoutfh = "QQFH";
	    }
	    
	    foreach my $readname ( @{$patientsort{$patient}}) {

		print $outfh $header{$readname}."\n";
		
		my @nextseq = split / */, $seq{$readname};
		$linecount = 0;

		$DEBUG && print LOG "DEBUG: length of sequence $readname is ",length($seq{$readname}), " or ", scalar(@nextseq)," - printing to $outfh\n";

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

		print $qualoutfh $header{$readname},"\n";
		print $qualoutfh $qual{$readname};
	    }
	}
    }

    close PFH;
    close QFH;	       

    close PQFH;
    close QQFH;
}

close LOG;

exit;
### ### ### ##### ###

# my @numbers = keys %header;
# my $entries = @numbers;

# for (my $i = 0; $i < $entries; $i++) {   
#     $pick_next = int(rand($entries - $i));

#     push @new_order, $numbers[$pick_next];
#     splice @numbers, $pick_next, 1;
# }

# if (@numbers > 0 ) {
#     print "Oops, ",scalar(@numbers)," numbers left in randomisation queue!";
# }



# my $setfh; 
# open PFH, ">$fastafile.p.set";
# open QFH, ">$fastafile.q.set";

# my $outfh = "PFH";
# my $pcount = 0;
# my $qcount = 0;

# while ( my $next = shift @new_order ) {

#     print $outfh $header{$next}, "\n";

#     my $linecount=0;
#     my @nextseq = split (/ */,$seq{$next});

#     $DEBUG && print "DEBUG: length of sequence $next is ",length($seq{$next}), " or ", scalar(@nextseq),"\n";       
    
#     while( $nextbase = shift @nextseq ) {
# 	$linecount++;
# 	print $outfh $nextbase;
# 	if($linecount == 60) {
# 	    print $outfh "\n";
# 	    $linecount = 0;
# 	}

#     }
#     if ($linecount != 0) {
# 	print $outfh "\n";
#     }

#     if($pcount/$entries > $set1_fract ) {
# 	$outfh = "QFH";
# 	$qcount++;
#     } else {
# 	$pcount++;
#     }
# }

# close PFH;
# close QFH;

# print STDERR "Order randomised. Wrote $pcount entries to $fastafile.p.set and $qcount to $fastafile.q.set.\n";
