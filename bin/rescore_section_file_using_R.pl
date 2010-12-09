#!/usr/bin/perl -w

my @clargs = @ARGV;
my $DEBUG = 0;
my $WARNING = 1;

my $SE_UKS_COMPAT = 1;

my $outfile = "";
my $infile = "";

my $count_one_dominant_per_isolate = 0;

sub usage;

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

if( ! -e $infile ) {
    $WARNING && print STDERR "FATAL: Could not open infile $infile.\n";
    usage();
}

$outfile = $infile.".tab";
if ( -e $outfile ) {
    $WARNING && print STDERR "Outfile $outfile exists. Cowardly refusing to overwrite.\n";
    usage();
}

open OUTTAB, ">$outfile";

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

#	    ($seed_id, $actual_members) = ( $row =~ /Seed ([\w\d]+_\d+_\d+) - (\d+) members/);
	    ($seed_id, $actual_members) = ( $row =~ /^Seed (\S+) - (\d+) members/);

	    # reset per isolate counter for each new motif (used by OPI model)
	    %counted_isolate = ();
	    $co_dominants_in_motif = 0;
	    $se_dominants_in_motif = 0;

#	    $motif_bulkout[$motif_nr] .= $row."\n"; # copy and paste
	    
	} elsif ( $row =~ m/(\s+)fraction/ ) { # SEfraction
	    # captures the "+" sense of the type_over-representation

	    my $over_type = $1;
	    if ($se_type ne $over_type) {
		# assume co_type equals over-type
		if( $co_type ne $over_type ) {
		    $WARNING && print "Found fraction over-rep type $over_type, but it is not equal to either $se_type or $co_type.\n";
		    # call as fatal?!
		} 

		$co_type = $se_type;
		$se_type = $over_type;
		
	    }
#	    my ($old_se_dominants_in_motif, $old_co_dominants_in_motif, $se_frac, $pval, $se_over) = ( $row =~/SE\s+(\d+)\s+CO\s+(\d+)\s+SEfraction ([.0-9]+) p ([.0-9]+) ([-+]{1})/ );

	    $DEBUG && print STDERR "OLD STATS ROW: $row\n";
	    # we can now compute and give a new se-skew and p-value for this motif

	    if($isolatetype_dominants_in_motif{$se_type} + $isolatetype_dominants_in_motif{$co_type} > 0 ) {
#		my $pval;
#		my $se_over;
#		$DEBUG && print STDERR "DEBUG: Calling dhyper for seed $seed_id ($actual_members actual members) with M=$se_dominants, N=$co_dominants, m=$se_dominants_in_motif, n=$co_dominants_in_motif.\n";
		$DEBUG && print STDERR "DEBUG: Calling dhyper for seed $seed_id ($actual_members actual members) with M=".($isolatetype_dominants{$se_type}).", N=".($isolatetype_dominants{$co_type}).", m=".($isolatetype_dominants_in_motif{$se_type}).", n=".($isolatetype_dominants_in_motif{$co_type})."\n";

		print OUTTAB "$seed_id\t$actual_members\t".($isolatetype_dominants{$se_type})."\t".($isolatetype_dominants{$co_type})."\t".($isolatetype_dominants_in_motif{$se_type})."\t".($isolatetype_dominants_in_motif{$co_type})."\n";

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

	    ($motif_id[$motif_nr],$motif_start,$motif_end, $fractcount) = ( $row =~ m/^\s*([\w\d]+)\s+(\d+)\s+(\d+)\s+\w+\s*(.*)\s*$/ );

	    $motif_start = $motif_start -1;
	    $motif_end = $motif_end -1;

	    if ($fractcount ne "") {
	    
		if($count_one_dominant_per_isolate == 1) {

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
			    $DEBUG && print STDERR "OPI model active: ignoring one $current_count.\n";
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
	    }
#	    $motif_bulkout[$motif_nr] .= $row."\n"; # copy and paste - no need to edit individual motif members; only final tally

	}
    } elsif ( $row =~ /From fraction-count: \w+ dominants: \d+ \w+ dominants: \d+\./)  {
	# also want to catch the totals
	($co_type, $co_dominants, $se_type, $se_dominants) = ( $row =~ /From fraction-count: (\w+) dominants: (\d+) (\w+) dominants: (\d+)\./);
    } elsif ($row =~ /^Found (\d+) motifs\./ ) {

	$DEBUG && print STDERR $row,"\n";

	($n_motifs) = ($row =~ /^Found (\d+) motifs\./);
	$motifs_ahead = 1;

	@motif_pval = 1 x $n_motifs;
#	@motif_bulkout = "" x $n_motifs;
	@motif_id = "" x $n_motifs;
    } else {
#	$motifs_only || print OUTTAB $row,"\n";
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

close SECTION;
close OUTTAB;

#
open RFILE, ">$outfile.r";
print RFILE "rm2<-read.table(\"$outfile\") ; \n";
print RFILE "sink(\"$outfile.hyp\") ; \n";
print RFILE 'for(i in (1:length(rm2$V5))) print(sprintf("%s %i %i %i %i %f", rm2$V1[i], rm2$V3[i], rm2$V4[i], rm2$V5[i],rm2$V6[i],phyper(m=rm2$V3[i], n=rm2$V4[i], q=ifelse(rm2$V3[i]/(rm2$V4[i]+rm2$V3[i])>rm2$V5[i]/(rm2$V5[i]+rm2$V6[i]), rm2$V5[i], rm2$V5[i]-1), k=rm2$V5[i]+rm2$V6[i], lower.tail= rm2$V3[i]/(rm2$V4[i]+rm2$V3[i])>rm2$V5[i]/(rm2$V5[i]+rm2$V6[i]) ))) ;',"\n";
print RFILE "sink()\n";
close RFILE;

system("R --no-save --slave < $outfile.r");

open HYP, "<$outfile.hyp";

while ($row = <HYP>) {

    $row =~ /\[\d+\]\s+\"(\S+)\s+\d+\s+\d+\s+\d+\s+\d+\s+([^\"]+)\"/;
    $seed = $1; 
    $hypergeoP = $2;

    $newP{$seed} = $hypergeoP;

}

close HYP;

open SECTION, "<$infile";

$seed_id = "";
$actual_members = "";
my $preserve_info = "";

$motif_nr = 0;
$motifs_ahead = 0;

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
		print STDERR "WARNING: bulk out contained $motif_bulkout[$motif_nr].\n";
		print STDERR "WARNING: Resolving this by deleting last motif.\n";

		$motif_bulkout[$motif_nr] = "";
	    }
	    
	    $last_motif_closed = 0;

#	    ($seed_id, $actual_members) = ( $row =~ /Seed ([\w\d]+_\d+_\d+) - (\d+) members/);
	    ($seed_id, $actual_members) = ( $row =~ /^Seed (\S+) - (\d+) members/);

	    # reset per isolate counter for each new motif (used by OPI model)
	    $motif_bulkout[$motif_nr] .= $row."\n"; # copy and paste
	} elsif ( $row =~ m/SEfraction/ ) {
	    #	    my ($old_se_dominants_in_motif, $old_co_dominants_in_motif, $se_frac, $pval, $se_over) = ( $row =~/SE\s+(\d+)\s+CO\s+(\d+)\s+SEfraction ([.0-9]+) p ([.0-9]+) ([-+]{1})/);
	    $DEBUG && print STDERR "OLD STATS ROW: $row\n";

	    $motif_bulkout[$motif_nr] .= $row." P ".$newP{$seed_id}."\n";

	    $motif_pval[$motif_nr] = $newP{$seed_id};

	    $last_motif_closed = 1;

	    $motif_nr++;
	} else {
	    #	    $DEBUG && print "GOT $row with motifs_ahead.\n";
	    
	    $motif_bulkout[$motif_nr] .= $row."\n"; # copy and paste - no need to edit individual motif members; only final tally
	    
	}
    } elsif ( $row =~ /From fraction-count: CO dominants: \d+ SE dominants: \d+\./)  {
	# also want to catch the totals
	($co_dominants, $se_dominants) = ( $row =~ /From fraction-count: CO dominants: (\d+) SE dominants: (\d+)\./);
	$preserve_info = $row."\n";
    } elsif ( $row =~ /^Find.+regions/ or $row =~ /^Lastly.+parts\./ or $row =~ /^Found.+regions/ or $row =~ /^Indicator peak/) {
	$preserve_info .= $row."\n";
    } elsif ( $row =~ /^Found (\d+) motifs\./ ) {

	$DEBUG && print STDERR $row,"\n";

	($n_motifs) = ($row =~ /^Found (\d+) motifs\./);
	$motifs_ahead = 1;

	@motif_pval = 1 x $n_motifs;
	@motif_bulkout = "" x $n_motifs;
	@motif_id = "" x $n_motifs;
    } else {

	$motifs_only || ($preserve_info .= $row."\n");
    }

}

($motif_nr == $n_motifs) || print STDERR "WARNING: found $motif_nr motifs but original file claimed $n_motifs motifs.\n";

open OUTFILE, ">$infile.newp";

print OUTFILE $preserve_info, "\n";
print OUTFILE "Found $n_motifs motifs.\n\n";

my @order = sort { $motif_pval[$a] <=> $motif_pval[$b] } (0..($motif_nr-1));

foreach my $i ( @order ) {
    print OUTFILE $motif_bulkout[$i];
    print OUTFILE "\n";
}

print OUTFILE "\nDone. ;-)\n";

sub usage {
    print "USAGE: rescore_section_file_using_R.pl\t<--in section>\n";
    print "\t\t\t\t[--countmodel oneperisolate]\t\tcount one dominant per isolate\n";
    print "\t\t\t\t[--motifs_only]\t\tRetain only the motif summary\n";
    exit(1);
}
