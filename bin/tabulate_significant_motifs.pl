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
open OUTTAB, ">$outfile";

my @motif_id;

my $motif_nr = 0;

my $motifs_ahead = 0;

my $last_motif_closed = 1;

my %motif_startcol;
my %motif_endcol;

# needs: section_file_name, count_one_dominant_per_isolate 
open SECTION, $infile;

my $seed_id;
my $actual_members;
my $seed_seq; 

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
	    $DEBUG && print STDERR "In $seed_id with $actual_members members.\n";
	    # reset per isolate counter for each new motif (used by OPI model)
#	    %counted_isolate = ();
#	    $co_dominants_in_motif = 0;
#	    $se_dominants_in_motif = 0;

#	    $motif_bulkout[$motif_nr] .= $row."\n"; # copy and paste
	    
	} elsif ( $row =~ m/SEfraction/ ) {
	    
	    my ($se_dominants_in_motif, $co_dominants_in_motif, $se_frac, $se_over, $pval) = ( $row =~/SE\s+(\d+)\s+CO\s+(\d+)\s+SEfraction ([.0-9]+) p [.0-9]+ ([-+]{1}) P ([.0-9Na]+)/ );

#	    $DEBUG && print STDERR "OLD STATS ROW: $row\n";
	    # we can now compute and give a new se-skew and p-value for this motif

#	    if($se_dominants_in_motif + $co_dominants_in_motif > 0 ) {
#		my $pval;
#		my $se_over;
#		$DEBUG && print STDERR "DEBUG: Calling dhyper for seed $seed_id ($actual_members actual members) with M=$se_dominants, N=$co_dominants, m=$se_dominants_in_motif, n=$co_dominants_in_motif.\n";

	    my ($startcol, $endcol);

	    if( defined($motif_startcol{$seed_id}) ) {
		$startcol = $motif_startcol{$seed_id};
		$endcol = $motif_endcol{$seed_id};
	    }

	    my $skew = ( $se_dominants_in_motif / $se_dominants - $co_dominants_in_motif / $co_dominants ) / ( $se_dominants_in_motif / $se_dominants + $co_dominants_in_motif / $co_dominants );


	    $DEBUG && print STDERR "$seed_seq\t$seed_id\t$actual_members\t$se_dominants_in_motif\t$co_dominants_in_motif\t$pval\t".sprintf("%.6f",$skew)."\t$startcol\t$endcol\n";

	    print OUTTAB "$seed_seq\t$seed_id\t$actual_members\t$se_dominants_in_motif\t$co_dominants_in_motif\t$pval\t".sprintf("%.6f",$skew);
	    if( defined($motif_startcol{$seed_id}) ) {
		print OUTTAB "\t$startcol\t$endcol";
	    }
	    print OUTTAB "\n";	    

	    $seed_seq ="";

#	($pval, $se_over) = dhyper($se_dominants, $co_dominants, $se_dominants_in_motif, $co_dominants_in_motif);
#		$motif_pval[$motif_nr] = $pval;
##	$motif_seover[$motif_nr] = $se_over;
#		$motif_bulkout[$motif_nr] .= sprintf "SE %3d CO %3d SEfraction %.3f p %.4f %s\n", $se_dominants_in_motif, $co_dominants_in_motif, ($se_dominants_in_motif / ( $se_dominants_in_motif+$co_dominants_in_motif)), $pval, $se_over;
#		$DEBUG && printf STDERR "SE %3d CO %3d SEfraction %.3f p %.4f %s\n", $se_dominants_in_motif, $co_dominants_in_motif, ($se_dominants_in_motif / ( $se_dominants_in_motif+$co_dominants_in_motif)), $pval, $se_over;
#	exec("echo |R")

	#    }
	    
	    $last_motif_closed = 1;

	    $motif_nr++;
	} else {
	    
#	    $DEBUG && print "GOT $row with motifs_ahead.\n";
	    
	    my $motif_start;
	    my $motif_end;
	    my $fractcount;
	    my $motif_member_sequence; 
	    ($motif_id[$motif_nr],$motif_start,$motif_end, $motif_member_sequence, $fractcount) = ( $row =~ m/^\s*([\w\d]+)\s+(\d+)\s+(\d+)\s+(\w+)\s*(.*)\s*$/ );
	    
	    $DEBUG && print STDERR "Motif member ",join(" ",$motif_id[$motif_nr],$motif_start,$motif_end, $motif_member_sequence, $fractcount),"\n";

	    my $full_id = $motif_id[$motif_nr]."_".$motif_start."_".$motif_end;
	    if ( $full_id eq $seed_id ) {
		$seed_seq = $motif_member_sequence;
	    }

	    $motif_start = $motif_start -1;
	    $motif_end = $motif_end -1;
	    	    
#	    $motif_bulkout[$motif_nr] .= $row."\n"; # copy and paste - no need to edit individual motif members; only final tally

	}
    } elsif ( $row =~ /From fraction-count: CO dominants: \d+ SE dominants: \d+\./)  {
	# also want to catch the totals
	($co_dominants, $se_dominants) = ( $row =~ /From fraction-count: CO dominants: (\d+) SE dominants: (\d+)\./);
	print OUTTAB "Total\tSE\t$se_dominants\tCO\t$co_dominants\n";
    } elsif ($row =~ /^Strict section from (\d+) to (\d+) for seed\s+([\w\d_]+)\./) {

	#
	# columns
	
	my $startcol = $1;
	my $endcol = $2;
	my $seedid = $3;
	
	$motif_startcol{$seedid} = $startcol;
	$motif_endcol{$seedid} = $endcol;

	$DEBUG && print "Motif $seedid in columns $startcol - $endcol.\n";

    } elsif ($row =~ /^Found (\d+) motifs\./ ) {

	$DEBUG && print STDERR $row,"\n";

	($n_motifs) = ($row =~ /^Found (\d+) motifs\./);
	$motifs_ahead = 1;

	@motif_id = "" x $n_motifs;
    } else {
#	$motifs_only || print OUTTAB $row,"\n";
    }
        
}

# first mark columns where id drops below treshold 

close SECTION;
close OUTTAB;

sub usage {
    print "USAGE: tabulate_significant_motifs.pl\t<--in section>\n";
    exit(1);
}
