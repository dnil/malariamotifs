#!/usr/bin/perl -w

use strict; 

my $DEBUG = 0;

my @primers = ( "AF", "NDBL", "NFBR");
my @stages = ("R","T");
my @seqprimers = ("f","r");

my $primerpattern = join('|', @primers);
my $stagepattern = join("", @stages);
my $seqprimerpattern = join("", @seqprimers);

$DEBUG && print "$primerpattern, $stagepattern, $seqprimerpattern\n";

my @contig; # array of contig-hashes
my @isolates; # book keep isolates for later ease of sorting

my $ignore_seqprimer = 0;

my $contig = "OOPS";
my %count;
my %seen_mate;

while(my $r= <STDIN>) {
    chomp $r;

    if( $r =~ /^CO\s+(\S+)\s+(\d+)\s+(\d+)\s+\d+\s+([UC])/) {
	# new contig
	$contig = $1;

#	$$contig[$currentcontignr]{name} = $1;
#	$$contig[$currentcontignr]{length} = $2;
#	$$contig[$currentcontignr]{n_reads} = $3;
#	$$contig[$currentcontignr]{direction} = $4;

    }

    if ( $r =~ /^AF\s+(\S+)/) {
	my $read = $1;
	
	$read =~ /^([$seqprimerpattern]{1})(\S+?)($primerpattern)([$stagepattern])\d*\.(\d+)\./;
	
	my $seqprimer = $1;
	my $isolate = $2;
	my $primer = $3;
	my $stage = $4;
	my $clone = $5;
	
#	$seqprimer && $isolate && $primer && $stage || 	print "$read\t$seqprimer - $isolate - $primer - $stage\n";
 	$DEBUG && print "$read\t$seqprimer - $isolate - $primer - $stage - $clone\n";

	if(!$ignore_seqprimer && defined($seen_mate{$isolate.$primer.$stage.$clone}) && $seen_mate{$isolate.$primer.$stage.$clone} == 1) {
	    # control for f-r pair
	} else {
	    $seen_mate{$isolate.$primer.$stage.$clone} = 1;

	    if( defined( ${$count{$contig}}{contig_total} ) ) { ${$count{$contig}}{contig_total}++ } else { ${$count{$contig}}{contig_total} = 1 };

	    # contig & primer
	    if( defined( ${${$count{$contig}}{$primer}}{$stage} ) ) { ${${$count{$contig}}{$primer}}{$stage}++ } else { ${${$count{$contig}}{$primer}}{$stage} = 1 };
	    if( defined( ${${$count{$contig}}{$primer}}{contigtotal_for_primer} ) ) { ${${$count{$contig}}{$primer}}{stagetotal_for_primer}++ } else { ${${$count{$contig}}{$primer}}{stagetotal_for_primer} = 1 };
	    
	    # contig & stage
	    if( defined( ${${$count{$contig}}{$stage}}{$isolate} ) ) { ${${$count{$contig}}{$stage}}{$isolate}++ } else { ${${$count{$contig}}{$stage}}{$isolate} = 1 };
	    if( defined( ${${$count{$contig}}{$stage}}{contigtotal_for_stage} ) ) { ${${$count{$contig}}{$stage}}{contigtotal_for_stage}++ } else { ${${$count{$contig}}{$stage}}{contigtotal_for_stage} = 1 };
#	    if( defined( ${$count{$contig}}{$stage}) ) { ${$count{$contig}}{$stage}++ } else { ${$count{$contig}}{$stage} = 1 };

#	    if ( defined( $count{isolate}{}{'stagetotal'}) ) { $count{stagetotal_for_primer}++ } else { $count{stagetotal_for_primer} = 1 };
		 
#	    if( defined( ${$count{$contig}}{primer_total} ) ) { ${$count{$contig}}{primer_total}++ } else { ${$count{$contig}}{primer_total} = 1 };



#	    if( defined( ${$count{$contig}}{stage_total} ) ) { ${$count{$contig}}{stage_total}++ } else { ${$count{$contig}}{stage_total} = 1 };

	    if( defined(${$count{$isolate}}{isolate_total}) ) { ${$count{$isolate}}{isolate_total}++ } else { ${$count{$isolate}}{isolate_total} = 1; push @isolates, $isolate;};
	    if( defined(${${$count{$isolate}}{$stage}}{isolatetotal_for_stage}) ) { ${${$count{$isolate}}{$stage}}{isolatetotal_for_stage}++ } else { ${${$count{$isolate}}{$stage}}{isolatetotal_for_stage} = 1 };
	}
    }    
}

foreach $contig (keys %count) {
    if($contig =~ /Contig/) {
	
	print "CONTIG $contig\t", ${$count{$contig}}{contig_total};
	
	# subtotals..
	foreach my $stage (@stages) {
	    if(defined(${${$count{$contig}}{$stage}}{contigtotal_for_stage})) { 
		print "\t$stage=",${${$count{$contig}}{$stage}}{contigtotal_for_stage};
	    } else {
		${${$count{$contig}}{$stage}}{contigtotal_for_stage} = 0;
	    }
	    

	    foreach my $isolate (keys %{ ${$count{$contig}}{$stage} }) {
		if($isolate ne "contigtotal_for_stage") {
		    print "\t$isolate.$stage=".${${$count{$contig}}{$stage}}{$isolate};		    
		}
	    }
	}
	
	foreach my $primer (@primers) {
	    if(defined(${${$count{$contig}}{$primer}}{contigtotal_for_primer})) {
		print "\t$primer=".${${$count{$contig}}{$primer}}{contigtotal_for_primer};
	    } else {
		${${$count{$contig}}{$primer}}{contigtotal_for_primer} = 0;
	    }
	}

	foreach my $stage (@stages) {
	    foreach my $primer (@primers) {
		if ( defined(${${$count{$contig}}{$primer}}{$stage}) ) {
		    print "\t".$stage."&".$primer,"=",${${$count{$contig}}{$primer}}{$stage};
		} else {
		    ${${$count{$contig}}{$primer}}{$stage} = 0;
		}
	    }	   
	}

    print "\n";
    } else {
	my $isolate = $contig;
	if(defined(${$count{$isolate}}{isolate_total})) {
	    print "ISOLATE $isolate ", ${$count{$isolate}}{isolate_total};
	} else {
	    ${$count{$isolate}}{isolate_total} = 0;
	}

	foreach my $stage (@stages) {
	    if (defined( ${${$count{$isolate}}{$stage}}{isolatetotal_for_stage} ) ) { 
		print "\t$stage=", ${${$count{$isolate}}{$stage}}{isolatetotal_for_stage};
	    } else {
		${${$count{$isolate}}{$stage}}{isolatetotal_for_stage} = 0;
	    }
	}
	print "\n";
    }
}

# gör också en lista sorterad på isolate, stage och sedan numeriskt på isolate.stage

my @contigs = keys %count; 

foreach my $isolate ( @isolates ) { # assuming all included isolates are available for the first stage
    $isolate eq "contigtotal_for_stage" && next;
    foreach my $stage (@stages) {
	
	# filter non-stage-isolate contigs (0 members)
	
	my @si_contigs = sort { ${${$count{$b}}{$stage}}{$isolate} <=> ${${$count{$a}}{$stage}}{$isolate} } grep defined(${${$count{$_}}{$stage}}{$isolate}), @contigs;

	foreach my $contig (@si_contigs) { 
	    if( !defined (${${$count{$contig}}{$stage}}{$isolate}) || (${${$count{$contig}}{$stage}}{$isolate} == 0)) {
		# next
	    } else {
		print "DOMINANCE $isolate\t$stage\t$contig\t",${${$count{$contig}}{$stage}}{$isolate},"\t",${${$count{$isolate}}{$stage}}{isolatetotal_for_stage},"\t",${$count{$contig}}{contig_total};
		foreach my $primer (@primers) {
		    defined(${${$count{$contig}}{$primer}}{$stage}) && 
			print "\t",$primer,"=",${${$count{$contig}}{$primer}}{$stage};		    
		}
		print "\n";
	    }

	}
    }
}


#ok, så, vad du också skulle behöva räkna ner är isolate.stage - totalen för varje contig. Skit i primer så länge.
# och så vill man veta hackordningen i varje isolate.stage - vilken contig har flest?

# USage: sed -e 's/D7AH1S2/AH1S2/g' stagecomp.061103.screen.plq.seq.ace |/home/daniel/malaria.0606/malariamotifs/bin/print_isolate_list_from_ace_file_stagecomp.pl  |sort -k1,1 -k2.8n|less
# You may also want to 
# |grep "^CONTIG\|^ISOLATE" |sort -k1,1 -k2.8n > stagecomp.061103.ace.summary_v2
# |grep "DOMINANCE" >> stagecomp.061103.ace.summary_v2
