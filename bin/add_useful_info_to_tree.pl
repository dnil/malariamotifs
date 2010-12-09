#!/usr/bin/perl -w

# count_fractions_in_ace_file.pl < se31reads.fasta.ace | sort -k2,2nr -k5,5nr

my $isolate_sort = 0;
my $DEBUG = 0;
my $WARNING = 1;

my $dominant_cutoff = 4;
my $dominance_exclusion_gbk_ids_file= "";
my $dominance_exclusion_contig_ids_file = "";
my $exclude = 0;
my $blasthitgff_file = "";

my $dominance_model = "sum"; # sum of primer contributions determines dominance order by default

# CLI

while (my $arg = shift @ARGV) {
    if($arg eq "--debug") {
	$DEBUG = 1;
    } elsif($arg eq "--typesort") {
	$typesort_file = shift @ARGV;	
    } elsif($arg eq "--fractioncount") {
	$fractioncount_file = shift @ARGV;
    } elsif($arg eq "--dominant_cutoff") {
	$dominant_cutoff = shift @ARGV;
    } elsif($arg eq "--blasthitgff") {
	$blasthitgff_file = shift @ARGV;
    } elsif($arg eq "--exclude_from_dominance") {
	$dominance_exclusion_gbk_ids_file = shift @ARGV;
    } elsif($arg eq "--exclude_do_contigs") {
	$dominance_exclusion_contig_ids_file = shift @ARGV;
    } elsif($arg eq "--dominance") {
	$dominance_model = shift @ARGV;
	if( $dominance_model eq "sum" or $dominance_model eq "sum_or_af_top") {
	    $DEBUG && print STDERR "DEBUG: Using dominance model $dominance.\n";	    
	} else {
	    print STDERR "WARNING: Unrecognised dominance model $dominance. Using sum model instead.\n";
	    $dominance_model = "sum";
	}       
    } elsif($arg eq "--tree") {
	$tree_file = shift @ARGV;
    } elsif ($arg eq "--no-warn") {	
	$WARNING = 0;
    } else {
	print STDERR "WARNING: Unrecognised command line parameter $arg. Continuing anyway.\n";
    }
}

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

open(TYPESORT, "<$typesort_file");
open(FRACTCOUNT, "<$fractioncount_file");
open(TREE, "<$tree_file");
open(BLASTGFF, "<$blasthitgff_file");

my %typesort;

while(<TYPESORT>) {
    chomp;
    my (@row) = split /\s+/,$_;

    if( exists($typesort{$row[0]}) ) {
	$typesort{$row[0]} .= "_".$row[1]."_".$row[2];
    } else {
	$typesort{$row[0]} = $row[1]."_".$row[2];
    }
}


my %blasthit;
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
    
    my $subjectname;
    if ($comment[0] =~ m/b\|([\w\d\.]+)\|/) {	
	$subjectname = $1;
    } else {
	($subjectname) = ($comment[0] =~ m/([\w\d\_\.]+)\|mRNA\|/);
    }

    my ($expect) = ($comment[1] =~ m/Expect = ([\deE\.\-]+)/);
    my ($querycoverage) = ($comment[1] =~ m/Querycoverage = ([\d\.]+)/);
    my ($idfreq) = ($comment[1] =~ m/Ids = \d+\/\d+ \(([\d\.]+)\)/);

    if(!defined($lowestexpect{$contigname}) || ($expect < $lowestexpect{$contigname})) {	
	$blasthit{$contigname} = "_".$subjectname."_$querystart+".$queryend."_$querycoverage"."_".$idfreq;
	$lowestexpect{$contigname} = $expect;
    }

    $DEBUG && print STDERR $contigname." ".$blasthit{$contigname}." ".$lowestexpect{$contigname}."\n";    
}

if($exclude == 1) {
    
    # check exclusion id list
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
	foreach my $exclude_contig (@exclude_contig) {
	    $exclude{$exclude_contig} = 1;
	    $DEBUG && print STDERR "Excluding $exclude_contig based on contig exclude list entry.\n";
	}
    }

}

my %isolate_fraction;
my %isolate_contig_primer_fraction;

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

my %fractcount;

foreach my $isolate (keys %isolate_fraction) {

    my (@contigs) = keys %${$isolate_fraction{$isolate}};

    if ( $exclude == 1 ) {
	my @include_contigs;
	foreach my $contig (@contigs) {
	    if ( defined( $exclude{$contig} ) ) {
		if( $exclude{$contig} == 1) {
		    
		    $DEBUG && print STDERR "Excluded contig $contig.\n";
		} elsif ($exclude{$contig} == 0) {
		    push @include_contigs, $contig;
		}
	    } else {
		# if contig didn't get a blast hit, mark it non-excluded =)
		$exclude{$contig} = 0;
		push @include_contigs, $contig;
	    }
	}
	@contigs = @include_contigs;
    }

    my $maxcount = ( @contigs < $dominant_cutoff ) ? @contigs  : $dominant_cutoff;

    $DEBUG && ($maxcount < $dominant_cutoff) && print STDERR "Got maxcount $maxcount from domcutoff $dominant_cutoff and ",scalar(@contigs)," contigs.\n";

    if( $dominance_model eq "sum" ) { 

	@contigs = sort { ${$isolate_fraction{$isolate}}->{$b} <=> ${$isolate_fraction{$isolate}}->{$a}} @contigs;

#    $DEBUG && print "DEBUG: $isolate maxcount $maxcount\n";

	for( my $i = 0; $i < $maxcount ; $i++) {

	    if( exists($fractcount{$contigs[$i]}) ) { 
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
	    foreach my $primer ( "af","nfbr","ndbl" ) {
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
    
    } elsif ( $dominance_model eq "primer_allaboard" ) {
	
	# need per isolate n contig fracs - present in frac_count file
	
	my %contigs_primer;

	my (@contigs) = keys %${$isolate_fraction{$isolate}};
	
	foreach my $contig ( @contigs ) {
	    foreach my $primer ( keys %{${$isolate_contig_primer_fraction{$isolate}}->{$contig}} ) {
		if ( ${${$isolate_contig_primer_fraction{$isolate}}->{$contig}}->{$primer} ) {
		    push @{$contigs_primer{$primer}}, $contig;
		}
	    }
	}

	foreach my $primer ( keys %contigs_primer ) {
	    
	    @{$contigs_primer{$primer}} = sort { ${${$isolate_contig_primer_fraction{$isolate}}->{$b}}->{$primer} <=> ${${$isolate_contig_primer_fraction{$isolate}}->{$a}}->{$primer} } @{$contigs_primer{$primer}};
	    
	}

	$dominants = 0;

	my $top = "";
	my $top_dominant_is_consistent = 0;
	foreach my $primer ( keys %contigs_primer ) {

	    if ($top eq "") {
		$top = ${$contigs_primer{$primer}}[0];
	    } elsif ($top eq ${$contigs_primer{$primer}}[0] ) {
		# two in agreement are enough, if there were only two primers used
		$top_dominant_is_consistent = 1;
	    } else {
		# one disagreeing is enough to use other model
		$top_dominant_is_consistent = 0;
		last;
	    }
	}

	if($top_dominant_is_consistent == 1) {	   
	    $fractcount{$top} .= "_".$isolate.".1_".${$isolate_fraction{$isolate}}->{$top};
	    $dominants = 1;
	}
	
    }
    
}

while(my $row = <TREE>) {
    chomp $row;
    if($row =~ m/(Contig\d+)/) {
	my $contig = $1;
	my $addin = "";

	if(exists($fractcount{$contig})) {
	    $addin = "\_".$fractcount{$contig};
	}

	$addin .= "\_".$typesort{$contig};

	if(defined($blasthit{$contig})) {
	    $addin .= "\_".$blasthit{$contig};
	} else {
	    $DEBUG && print STDERR "No hits for $contig.\n";
	}

	if( defined($exclude{$contig}) && ($exclude{$contig} == 1) ) {
	    $addin .= "\_"."dox"; 
	}
	#$row =~ s/$contig/$contig\_$fractcount{$contig}\_$typesort{$contig}/;


	$row =~ s/$contig/$contig$addin/;
    }
    
    print $row."\n"; 
}
 
   
