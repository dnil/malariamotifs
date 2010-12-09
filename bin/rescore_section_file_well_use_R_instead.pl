#!/usr/bin/perl -w

my @clargs = @ARGV;
my $DEBUG = 0;
my $WARNING = 1;

my $outfile = "";
my $infile = "";

my $count_one_dominant_per_isolate = 0;

sub dhyper;
sub over;
sub fac;
sub pfac;

my $motifs_only = 0;

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
    } elsif ($arg eq "--countmodel") {
	$arg = shift @clargs;
	if ($arg eq "oneperisolate") {
	    $count_one_dominant_per_isolate = 1;
	} else {
	    print "WARNING: unknown argument given to --countmodel.\n";	    
	}
    } elsif($arg eq "--motifs_only") {
	$motifs_only = 1;
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

$DEBUG && $motifs_only && print STDERR "DEBUG: printing mode set to motifs_only.\n";

my @motifs;
my @motif_pval;
my @motif_bulkout ;
my @motif_id;

my $motif_nr = 0;

my $motifs_ahead = 0;

my %counted_isolate = ();

my $last_motif_closed = 1;

# needs: section_file_name, count_one_dominant_per_isolate 
open SECTION, $infile;

my $seed_id;
my $actual_members;

while(my $row = <SECTION>) {

    chomp $row;

    if ($motifs_ahead) {

	if ( $row eq "" ) {
	    next;
	}

	if($row =~ /Done./ ) {
	    last;
	}

	if( $row =~ /^Seed/ ) {

	    if($last_motif_closed == 0) {
		print STDERR "WARNING: last motif did not terminate properly before next Seed entry. Please re-check motif finding & input data.\n";
#		print STDERR "WARNING: bulk out contained $motif_bulkout[$motif_nr].\n";
		print STDERR "WARNINIG Resolving this by deleting last motif.\n";

#		$motif_bulkout[$motif_nr] = "";
	    }
	    
	    $last_motif_closed = 0;

	    ($seed_id, $actual_members) = ( $row =~ /Seed (Contig\d+_\d+_\d+) - (\d+) members/);

	    # reset per isolate counter for each new motif (used by OPI model)
	    %counted_isolate = ();
	    $co_dominants_in_motif = 0;
	    $se_dominants_in_motif = 0;

#	    $motif_bulkout[$motif_nr] .= $row."\n"; # copy and paste
	    
	} elsif ( $row =~ m/SEfraction/ ) {
	    
#	    my ($old_se_dominants_in_motif, $old_co_dominants_in_motif, $se_frac, $pval, $se_over) = ( $row =~/SE\s+(\d+)\s+CO\s+(\d+)\s+SEfraction ([.0-9]+) p ([.0-9]+) ([-+]{1})/ );

	    $DEBUG && print STDERR "OLD STATS ROW: $row\n";
	    # we can now compute and give a new se-skew and p-value for this motif

	    if($se_dominants_in_motif + $co_dominants_in_motif > 0 ) {
#		my $pval;
#		my $se_over;
#		$DEBUG && print STDERR "DEBUG: Calling dhyper for seed $seed_id ($actual_members actual members) with M=$se_dominants, N=$co_dominants, m=$se_dominants_in_motif, n=$co_dominants_in_motif.\n";

		print "$seed_id\t$actual_members\t$se_dominants\t$co_dominants\t$se_dominants_in_motif\t$co_dominants_in_motif\n";

#	($pval, $se_over) = dhyper($se_dominants, $co_dominants, $se_dominants_in_motif, $co_dominants_in_motif);
#		$motif_pval[$motif_nr] = $pval;
##	$motif_seover[$motif_nr] = $se_over;
#		$motif_bulkout[$motif_nr] .= sprintf "SE %3d CO %3d SEfraction %.3f p %.4f %s\n", $se_dominants_in_motif, $co_dominants_in_motif, ($se_dominants_in_motif / ( $se_dominants_in_motif+$co_dominants_in_motif)), $pval, $se_over;
#		$DEBUG && printf STDERR "SE %3d CO %3d SEfraction %.3f p %.4f %s\n", $se_dominants_in_motif, $co_dominants_in_motif, ($se_dominants_in_motif / ( $se_dominants_in_motif+$co_dominants_in_motif)), $pval, $se_over;
#	exec("echo |R")

	    }
	    
	    $last_motif_closed = 1;

	    $motif_nr++;
	} else {
	    
#	    $DEBUG && print "GOT $row with motifs_ahead.\n";
	    
	    my $motif_start;
	    my $motif_end;
	    my $fractcount;

	    ($motif_id[$motif_nr],$motif_start,$motif_end, $fractcount) = ( $row =~ m/(Contig\d+)\s+(\d+)\s+(\d+)\s+\w+\s*(.*)\s*$/ );

	    $motif_start = $motif_start -1;
	    $motif_end = $motif_end -1;

	    if ($fractcount ne "") {
	    
		if($count_one_dominant_per_isolate == 1) {
		    while ( $fractcount =~ m/(co\d+)/g ) {
			$current_count = $1;
			if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {
			    $DEBUG && print STDERR "OPI model active: ignoring one $current_count.\n";
			} else {
			    $counted_isolate{$current_count} = 1;
			    $co_dominants_in_motif++;
			}
		    }
		
		    while ( $fractcount =~ m/(cp\d+)/g ) {			
			$current_count = $1;

			$current_count =~ s/cp/co/;
			if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {			    
			    $DEBUG && print STDERR "OPI model active: ignoring one $current_count.\n";
			} else {
			    $counted_isolate{$current_count} = 1;
			    $co_dominants_in_motif++;
			}
		    }
		    
		    while ( $fractcount =~ m/(uks\d+)/g ) {
			$current_count = $1;
			
			if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {			    
			    $DEBUG && print STDERR "OPI model active: ignoring one $current_count.\n";
			} else {
			    $counted_isolate{$current_count} = 1;
			    $se_dominants_in_motif++;
			}
		    }
		    
		    while ( $fractcount =~ m/(se\d+)/g ) {
			$current_count = $1;
			if(defined($counted_isolate{$current_count}) && $counted_isolate{$current_count} == 1) {			
			    $DEBUG && print STDERR "OPI model active: ignoring one $current_count.\n";
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
	    }

#	    $motif_bulkout[$motif_nr] .= $row."\n"; # copy and paste - no need to edit individual motif members; only final tally

	}
    } elsif ( $row =~ /From fraction-count: CO dominants: \d+ SE dominants: \d+\./)  {
	# also want to catch the totals
	($co_dominants, $se_dominants) = ( $row =~ /From fraction-count: CO dominants: (\d+) SE dominants: (\d+)\./);
    } elsif ( $row =~ /^Found (\d+) motifs\./ ) {

	print STDERR $row,"\n";

	($n_motifs) = ($row =~ /^Found (\d+) motifs\./);
	$motifs_ahead = 1;

	@motif_pval = 1 x $n_motifs;
#	@motif_bulkout = "" x $n_motifs;
	@motif_id = "" x $n_motifs;
    } else {
	$motifs_only || print $row,"\n";
    }
        
}

#($motif_nr == $n_motifs) || print STDERR "WARNING: found $motif_nr motifs but original file claimed $n_motifs motifs.\n";

#my @order = sort { $motif_pval[$a] <=> $motif_pval[$b] } (0..($motif_nr-1));

#foreach my $i ( @order ) {
#    print $motif_bulkout[$i];
#    print "\n";
#}

#print "\nDone. ;-)\n";

# first mark columns where id drops below treshold 

sub usage {
    print "USAGE: section_alignment.pl\t<--in section>\n";
    print "\t\t\t\t[--countmodel oneperisolate]\t\tcount one dominant per isolate\n";
    print "\t\t\t\t[--motifs_only]\t\tRetain only the motif summary\n";
    exit(1);
}
