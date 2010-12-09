#!/usr/bin/perl -w
#
# Parse ncbi-blast results, in particular MegaBlast hits

sub print_hit;

my $maxeval = -1;
my $DEBUG = 0;
my $WARNING = 1;

my $subjectHitBegin;
my $subjectHitEnd;
my $subjectHitStrand;
my $subjectName;
my $subjectLength; # currently unused..
my $queryHitBegin;
my $queryHitEnd;
my $queryHitStrand;
my $queryLength;
my $currentQuery;
my $queryName;
my $identities;
my $alignLength;
my $score;
my $expect;

my $coverpart = 0;
my $idst = 0;

my $lastQueryName = "first";
my $gffoutput = 1;

my $skim;

my $subjectref = 0;

while (my $arg = shift @ARGV) {
    if($arg eq "--eval" || $arg eq "-e") {
	$maxeval = shift @ARGV;
	if($maxeval eq "") {
	    print("Error! No max E value given.");
	    exit;
	}
    } elsif($arg eq "--no-gff" || $arg eq "-g") {
	$gffoutput = 0;
    } elsif($arg eq "--ids" || $arg eq "-i") {
	$idst = shift @ARGV;
	if($idst eq "") {
	    print("Error! No identity fraction value given.");
	    exit;
	} elsif ($idst > 1 || $idst< 0) {
	    print("Invalid identity fraction value ($coverpart) given.");
	    exit;
	}
    } elsif($arg eq "--cover" || $arg eq "-c") {
	$coverpart = shift @ARGV;
	if($coverpart eq "") {
	    print("Error! No coverage fraction value given.");
	    exit;
	} elsif ($coverpart > 1 || $coverpart < 0) {
	    print("Invalid coverage fraction value ($coverpart) given.");
	    exit;
	}
    } elsif($arg eq "--skim") {
	$skim = shift @ARGV;
	if($skim eq "") {
	    print("Error! No skim partial value given.");
	    exit;
	} elsif ($skim > 1 || $skim < 0) {
	    print("Invalid skim fraction value ($skim) given.");
	    exit;
	}		
    } elsif($arg eq "--subject") {
	$subjectref = 1;
    } elsif($arg eq "--debug") {
	$DEBUG = 1;
    } elsif($arg eq "--help" || $arg eq "-h" || $arg eq "-?") {
	usage();       
    }
}

sub usage {
    print STDERR "Usage: parselargemegablast.pl < blast_output > filtered_blast_hits_output\n";
    print STDERR "\t\t\t[--eval|-e max_e_value]\t\tFilter blast hits according to E-value.\n";
    print STDERR "\t\t\t[--cover|-c cover_part]\t\tFilter blast hits according to query coverage of hit, 0 =< cover_part <= 1\n";
    print STDERR "\t\t\t[--ids|-i cover_part]\t\tFilter blast hits according to identities of hit, 0 =< ids <= 1\n";
    print STDERR "\t\t\t[--skim skim_fraction]\t\tOutput randomly selected hits. \n\t\t\t\t\t\t\tThe percentage of hits returned equals skim_fraction 0<=skim_fraction<=1.\n";
    print STDERR "\t\t\t[--no-gff|-g]\t\tOutput in a simple non-gff format.\n";
    print STDERR "\t\t\t[--debug]\t\tEnable debug output.\n";
    print STDERR "\t\t\t[--help|-h|-?]\t\tDisplay this text.\n";
    exit;
}

# There is probably som better way to do this, but, hey, it works! =O)
my ($program);

# This version ok for multiple queries in one file.
# Ok for multi line subject-names.

my $nHits = 0;
my $nHitsCurrentQuery = 0;
my $filtered = 0;
my $inHit = 0;
my $inHitName = 0;
my $inQueryName = 0;
my $inDetail = 0;
my $detailField = 0;

my $subjectNamePart = "";
my $queryNamePart = "";

my $next_row_is_alignment = 0;

my $skipCurrentQuery = 0;

# nHits policy
# if you find a new hit, increase hit counter
# if old hit data is going away, print last hit.

while(<STDIN>) {
#    $ncbi_reply .= $_;

    s/\<a name \= \d+\>/\>/g; # Introduce > before hits, by substitution with the click-map name tag..
    s/\<.+?\>//g;

    if(/^MEGABLAST/) {
	$program = 'megablast';	
	$DEBUG && print STDERR "\nParsing (ncbi) $program output..\n";
	next;
    } elsif(/^TBLASTN/) {
	$program = 'tblastn';
	$DEBUG && print STDERR "\nParsing (ncbi) $program output..\n";
	next;
    } elsif(/^BLASTX/) {
	$program = 'blastx';
	$DEBUG && print STDERR "\nParsing (ncbi) $program output..\n";
	next;
    } elsif(/^BLASTN/) {
	$program = 'blastn';
	$DEBUG && print STDERR "\nParsing (ncbi) $program output..\n";
	next;
    } elsif(/^BLASTP/) {
	# untested & unimplemented
	$program = 'blastp';
	$DEBUG && print STDERR "\nParsing (ncbi) $program output..\n";
	next;
    }

    # TAKE CARE WHEN EDITING ANY PATTERN - A FEW OF THEM OCCUR IN SEVERAL PLACES, SO CHANGE THEM ALL!
    
    # Either we are looking at a detailed view of a hit, or we are scanning the hit-header.
    # The rest of the output is kindly enough different.. =)

    # how about printing query name, when two entries in a row with no hits were found and the output is not GFF? Care for it!
    # AND keep in mind that the last hit of the "previous" query may not have been printed.

    if($inQueryName) {
	if( /\(\s*[\d,]+\s+letters\)/ ) {
	    $inQueryName = 0;
	    ($queryLength) = /\(\s*([\d,]+)\s+letters\)/;
	    $queryLength=~s/,//;

	    if($queryLength == 0) {
		$WARNING && print STDERR "WARNING: ignoring any hits to qurrent query ($currentQuery) of length 0.\n"; 
		$skipCurrentQuery = 1;
	    }
	} else {
	    ($queryNamePart) = /^\s*(.+?)\s*$/; # remove initial/trailing ws
	    $currentQuery .= " $queryNamePart";
	}
	next;
    }

    if($inHit) {
	if($inHitName) {
	    if( /Length\s+\=\s+\d+/ ) {
		$inHitName = 0;
		($subjectLength) = /Length\s+\=\s+(\d+)/
	    } else {
		($subjectNamePart) = /^\s*(.+?)\s*$/; # remove initial/trailing ws
		$subjectName  .= " $subjectNamePart";
	    }
	}
	if($inQueryName) {
	    if( /\(\s*[\d,]+\s+letters\)/ ) {
		$inQueryName = 0;
	    } else {
		($queryNamePart) = /^\s*(.+?)\s*$/; # remove initial/trailing ws
		$currentQuery .= " $queryNamePart";
	    }
	}
	if($inDetail) {
	    if($next_row_is_alignment) {
		$alignment .= "$_\n";
		$next_row_is_alignment = 0;
	    }
	    if(/^Sbjct\:/) {
		if($detailField==0) {
		    # Get both hitBegin and hitEnd on first encounter, then only the hitEnds..
		    ($subjectHitBegin,$subjectHitEnd)=/^Sbjct\:\s+(\d+)\s*[\w*-]+\s+(\d+)/;
		} else { 
		    ($subjectHitEnd)=/^Sbjct\:\s+\d+\s*[\w*-]+\s+(\d+)/;
		}
		$alignment .= "$_\n";
	    } elsif(/^Query\:/) {
		# A new block of hit alignment rows was found.
		$detailField++;
		($queryHitEnd)=/^Query\:\s+\d+\s*[\w*-]+\s+(\d+)/;
		$alignment .= "\n$_\n"; #the extra newline only inside alignments, not before!
		$next_row_is_alignment = 1;
	    } elsif(/Score\s{1}\=/) {
		# Parse of alignment rows found a new hit on the same subject sequence.
		# remains unchanged
#		$subjectName[$nHits]=$subjectName;  

		print_hit;
		# found a new hit - increase counter
		$nHits++;
		$nHitsCurrentQuery++;

		$inDetail=0;
		
		($score)=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
#	  ($p)=/P\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
		($expect)=/Expect\(*\d*\)*\s*\=\s+([0-9\.e\-\+]+)/;
	    } elsif(/^\>/) {
		print_hit;
		$nHits++;
		$nHitsCurrentQuery++;

		# Hits supposedly begin with a row "> FASTA_NAME"
		($subjectName) = /^\>(.+?)\s*$/;

		# Parse of alignment rows found hit on naew subject sequence.
		$inDetail = 0;
		$inHitName = 1;

	    } elsif(/^Query\=/) {
		$skipCurrentQuery = 0; # initially assume querylen > 0, until checked

		# End of this query sequence. End of detail, end of hit.

		# print the last hit, since this is a new query sequence
		print_hit;
		$nHitsCurrentQuery = 0;

		$inDetail = 0;
		$inHit = 0;

		($currentQuery) = /^Query\=\s*(.+?)\s*$/;
		
		$inQueryName = 1;
	    } 
	    # in detail ends
	    #Parse a hit header..
	} elsif(/Score\s+\=/) {
	    ($score)=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
	    # Syntax of the P value varies btw runs in blastn.. *sigh*
	    ($expect)=/Expect\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
	} elsif (/Identities/) {
	    ($identities,$alignLength)=/Identities\s+\=\s+(\d+)\/(\d+)/;
	} elsif (/Frame\s+\=/) {
	    # strand is used in blastn & megablast parses
	    $queryHitStrand = '+';
	    $subjectHitStrand = /Frame\s+\=\s+([\d+-]+)/;
	} elsif (/Strand\s+\=/) {
	    # strand is used in blastn & megablast parses
	    ($queryHitStrand,$subjectHitStrand)=/Strand\s+\=\s+(\w+)\s+\/\s+(\w+)/;
	} elsif(/^Query\:/) {
	    $inDetail=1;
	    $detailField=0;
	    # If this is a gapped alignment, the aligned sequences may contain dashes for gaps.. 
	    ($queryHitBegin,$queryHitEnd)=/^Query\:\s+(\d+)\s*[\w*-]+\s+(\d+)/;
	    # Get both hitBegin and hitEnd on first encounter, later only the hitEnds..
	    $alignment .= "$_\n";
	    $next_row_is_alignment = 1;
	} 
    } elsif(/^Query\=/) { # notice: megablast/blastall *WILL* contain multiple queries..
	$skipCurrentQuery = 0; # initially assume querylen > 0, until checked
	
	if($lastQueryName ne 'first') { # if this isn't the first query
	    # print the last hit of last query since this is a new query sequence

	    print_hit;
	    $nHitsCurrentQuery = 0;

	    $inDetail = 0;
	    $inHit = 0;
	}

	($currentQuery) = /^Query\=\s*(.+?)\s*$/;
	$inQueryName = 1;
	$DEBUG && print STDERR "DEBUG: Query $currentQuery.\n";
    } elsif($skipCurrentQuery) { # the order is of essence: check this after Query= (although normally we'll find program name first), but certainly before checking for hits (>.+).
	next;
    } elsif(/^\>/) {

	# print previous hit

	# Apparently, the ORF-finding thing needs smallest-first type of coordinates?!
	print_hit;

	# found hit -- increase hit counter
	$nHits++;
	$nHitsCurrentQuery++;

	# Hits supposedly begin with a row "> FASTA_NAME"
	($subjectName)=/^\>(.+?)\s*$/;
	$inHitName = 1;

	$inHit=1;
	$alignment="";
    }
}

# print the last hit (and/or queryName)
print_hit;

if($nHits == 0) {
    print STDERR "WARNING: no hits found when parsing $program results.\n";
} else { 
    print STDERR "Found $nHits hits when parsing $program results.\n";
    $filtered && print STDERR "Filtered out $filtered of those according to criteria.\n";
}

# $n_hits = @{$hit{subjectName}};
# print "Found $n_hits hits.\n";

# init

# for($i=0;$i<$n_hits;$i++) {    
# }

sub print_hit {
    $DEBUG && print STDERR "$lastQueryName\n";

    # First, save the name of the new hit (not the one being printed)
    $queryName = $currentQuery;

    if ($lastQueryName ne $queryName) {
 	# output at first hit of this query (if query name has changed since last printed)
 	$lastQueryName = $queryName;
 	($gffoutput == 0) && print ">$lastQueryName\n";
     }
    
    if ($nHitsCurrentQuery == 0) {
	# print only query name, no hit
	return;
    }

    # Apparently, the ORF-finding thing needs smallest-first type of coordinates?!
    if($subjectHitBegin > $subjectHitEnd) {
	$tmp = $subjectHitBegin;
	$subjectHitBegin = $subjectHitEnd;
	$subjectHitEnd = $tmp ;
    }
    if($queryHitBegin > $queryHitEnd) {
	$tmp = $queryHitBegin;
	$queryHitBegin = $queryHitEnd;
	$queryHitEnd = $tmp ;
    }

    if($skim) {
	$randomint = int(rand(100));
	if( $randomint > $skim ) {
	    return; # discard due to skim.	    
	}
    }

    if($expect =~ /^e/) {
	# silly NCBI format output error
	($exp) = $expect =~ /^(e.+)\s*$/;
	$expect = "1$exp";
    }

# if subjectHits are needed:
#    if( $program eq 'tblastn' ) { # query length in nt, begin & end in aa coords.
#	$hitcover = sprintf("%.2f",($queryHitEnd - $queryHitBegin + 1) / $queryLength);	
#    } else {
#    }

    $hitcover = sprintf("%.2f",($queryHitEnd - $queryHitBegin + 1) / $queryLength);

    $idfrac = sprintf("%.2f", $identities/$alignLength);

    if ( ($maxeval > -1 && $expect > $maxeval) || ($coverpart > 0 && $hitcover < $coverpart) || ($idst > 0 && $idfrac < $idst ) ) {	
	# don't print sorely scoring hits.
	$filtered++;
    } else {
	if($program eq 'megablast' || $program eq 'blastn') {
	    $subjectStrand = "+";
	    if($subjectHitStrand eq "Minus") {
		$subjectStrand = "-";
	    }
	} elsif ($program eq 'tblastn' or $program eq 'blastx') {
	    $subjectStrand = "+";
	    if($subjectStrand eq '-1' || $subjectStrand eq '-2' || $subjectStrand eq '-3') {
		$subjectStrand = "-";
	    }
	} elsif ($program eq 'blastp') {
	    $subjectStrand = "+";
	    $subjectHitStrand = "+";
	    $queryHitStrand = "+";
	}

	if($gffoutput) {
	    ($gffQueryName) = ($lastQueryName =~ /^(.+)\s*/g);
	    if($subjectref == 0) {

		print "$gffQueryName\t$program\t".$program."hit\t$queryHitBegin\t$queryHitEnd\t.\t$subjectStrand\t.\tname $subjectName; comment Name = $subjectName\t, Score = $score, Expect = $expect, Ids = ",$identities,"/",$alignLength," (",$idfrac*100,"), Querycoverage = ",$hitcover*100,", Query hit $queryHitBegin-$queryHitEnd ($queryHitStrand), Subject hit $subjectHitBegin-$subjectHitEnd ($subjectHitStrand);\n";
	    } else {
		($subjectQueryName) = ($subjectName =~ /^(\S+)/g);
		# subject as reference!
		print "$subjectQueryName\t$program\t".$program."hit\t$subjectHitBegin\t$subjectHitEnd\t.\t$subjectStrand\t.\tname $subjectName; comment Name = $subjectName\t, Score = $score, Expect = $expect, Ids = ",$identities,"/",$alignLength," (",$idfrac*100,"), Querycoverage = ",$hitcover*100,", Query $gffQueryName hit $queryHitBegin-$queryHitEnd ($queryHitStrand), Subject hit $subjectHitBegin-$subjectHitEnd ($subjectHitStrand);\n";
	    }
	} else {
	    print "Name = $subjectName\t, Score = $score, Expect = $expect, Ids = ",$identities,"/",$alignLength," (",$idfrac*100,"%), Query hit $queryHitBegin-$queryHitEnd ($queryHitStrand), Coverage $hitcover, Subject hit $subjectHitBegin-$subjectHitEnd ($subjectHitStrand)\n";
	}
    }
    # NOTE: if reverse mapping to original ESTs is required, the trashed ESTs
    # must be accounted for!
}


#sub parse_megablast_reply {
#my $hit=shift;
#    print STDERR "\nParsing ncbi_reply..\n";
#    my $ncbi_reply=shift;

    # If HTML-formatted we translate the name tag into a fasta-header, and remove remaining tags
#    $DEBUG && print "DEBUG: $$ncbi_reply\n";

#    foreach $_ (split(/\n{1}/,$$ncbi_reply)) {

#      For($i=0; $i<$nHits;$i++) {
#  	$hit{subjectName}->[$i] = $subjectName[$i];
#  	$hit{subjectHitBegin}->[$i] = $subjectHitBegin[$i];
#  	$hit{subjectHitEnd}->[$i] = $subjectHitEnd[$i];
#  	$hit{subjectHitStrand}->[$i] = $subjectHitStrand[$i];
#  	$hit{queryName}->[$i] = $queryName[$i];
#  	$hit{queryHitBegin}->[$i] = $queryHitBegin[$i];
#  	$hit{queryHitEnd}->[$i] = $queryHitEnd[$i];
#  	$hit{queryHitStrand}->[$i] = $queryHitStrand[$i];
#  	$hit{identities}->[$i] = $identities[$i];
#  	$hit{alignLength}->[$i] = $alignLength[$i];
#  	$hit{score}->[$i] = $score[$i];
#  	$hit{expect}->[$i] = $expect[$i];
	
#  	$DEBUG && print "DEBUG: $queryName[$i] $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitStrand[$i]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ($subjectHitStrand[$i]) ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]\n";  
#      }

# }

#  sub parse_blastn_reply {
#      my $hit=shift;
#      my $ncbi_reply=shift;

#      # head
#      # ref, about query, 
#      # areamap
#      # Sequences producing significant alignments
#      # <PRE>
#      #  <a name = > </a><a href>id</a>name
#      # Length
#      # score, ids, strand  
#      # </PRE>
#      #<PRE> 
#      # Database:
    
#      # Ok, for a crude first attempt: remove HTML-tags. We lose information this way, 
#      # but on the other hand it is possible to use previously written WU-blast parsing code..

#      $$ncbi_reply=~s/\<a name \= \d+\>/\>/g; # Introduce > before hits, by substitution with the click-map name tag..
#      $$ncbi_reply=~s/\<.+?\>//g;
#      $DEBUG && print "DEBUG: $$ncbi_reply\n";

#      # IMPROVEMENT: Outdated parser version. Please introduce available code for long hit-names etc.

#      my $nHits=0;
#      my $inHit=0;
#      my $inDetail=0;
#      my $detailField=0;

#      my @subjectHitBegin;
#      my @subjectHitEnd;
#      my @subjectHitStrand;
#      my @subjectName;
#      my @queryHitBegin;
#      my @queryHitEnd;
#      my @queryHitStrand;
#      my $queryName;
#      my @identities;
#      my @alignLength;
#      my @score;
#      my @expect;
    
#      my $next_row_is_alignment = 0;

#      foreach $_ (split(/\n{1}/,$$ncbi_reply)) {    

#  	# TAKE CARE WHEN EDITING ANY PATTERN - A FEW OF THEM OCCUR IN SEVERAL PLACES, SO CHANGE THEM ALL!
	
#  	# Either we are looking at a detailed view of a hit, or we are scanning the hit-header.
#  	# The rest of the output is kindly enough different.. =)
#  	if($inHit) {
#  	    if($inDetail) { 
#  		if($next_row_is_alignment) {
#  		    $alignment[$nHits-1] .= "$_\n";
#  		    $next_row_is_alignment = 0;
#  		}
#  		if(/^Sbjct\:/) {
#  		    if($detailField==0) {
#  			# Get both hitBegin and hitEnd on first encounter, then only the hitEnds..
#  			($subjectHitBegin[$nHits-1],$subjectHitEnd[$nHits-1])=/^Sbjct\:\s+(\d+)\s+[atgcnxATGCNX-]+\s+(\d+)/;
#  		    } else { 
#  			($subjectHitEnd[$nHits-1])=/^Sbjct\:\s+\d+\s+[atgcnxATGCNX-]+\s+(\d+)/;
#  		    }
#  		    $alignment[$nHits-1] .= "$_\n";
#  		} elsif(/^Query\:/) {
#  		    # A new block of hit alignment rows was found.
#  		    $detailField++;
#  		    ($queryHitEnd[$nHits-1])=/^Query\:\s+\d+\s+[atgcnxATGCNX-]+\s+(\d+)/;
#  		    $alignment[$nHits-1] .= "\n$_\n"; #the extra newline only inside alignments, not before!
#  		    $next_row_is_alignment = 1;
#  		} elsif(/Score\s{1}\=/) {
#  		    # Parse of alignment rows found a new hit on the same subject sequence.
#  		    $subjectName[$nHits]=$subjectName[$nHits-1];
#  		    $inDetail=0;
#  		    $nHits++;
#  		    ($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
#  #	  ($p[$nHits-1])=/P\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
#  		    ($expect[$nHits-1])=/Expect\s*\=\s+([0-9\.e\-\+]+)/;
#  		} elsif(/^\>/) {
#  		    # Parse of alignment rows found hit on new subject sequence.
#  		    $inDetail=0;
#  		    # Hits supposedly begin with a row "> FASTA_NAME"
#  		    ($subjectName[$nHits])=/^\>(.+)/;
#  		    $nHits++;	     
#  		} 
#  		# in detail ends
#  		#Parse a hit header..     
#  	    } elsif(/Score\s{1}\=/) {
#  		($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
#  		# Syntax of the P value varies btw runs in blastn.. *sigh*
#  		($expect[$nHits-1])=/Expect\s+=\s+([0-9\.e\-\+]+)/;
#  	    } elsif (/Identities/) {
#  		($identities[$nHits-1],$alignLength[$nHits-1])=/Identities\s+\=\s+(\d+)\/(\d+)/;
#  	    } elsif (/Strand\s+\=/) {
#  		($queryHitStrand[$nHits-1],$subjectHitStrand[$nHits-1])=/Strand\s+\=\s+(\w+)\s+\/\s+(\w+)/;
#  	    } elsif(/^Query\:/) {
#  		$inDetail=1;
#  		$detailField=0;
#  		# If this is a gapped alignment, the aligned sequences may contain dashes for gaps.. 
#  		($queryHitBegin[$nHits-1],$queryHitEnd[$nHits-1])=/^Query\:\s+(\d+)\s+[atgcnxATGCNX-]+\s+(\d+)/;
#  		# Get both hitBegin and hitEnd on first encounter, later only the hitEnds..
#  		$alignment[$nHits-1] .= "$_\n";
#  		$next_row_is_alignment = 1;
#  	    }   
#  	} elsif(/^\>/) {
#  	    # Hits supposedly begin with a row "> FASTA_NAME"
#  	    ($subjectName[$nHits])=/^\>(.+)/;
#  	    $nHits++;
#  	    $inHit=1;
#  	    $alignment[$nHits-1]="";
#  	} elsif(/^Query\=/) {
#  	    # Actually just evaluated once to get query name..
#  	    ($queryName)=/^Query\=\s*(.+)/; 
#  #      print "Query $queryName.\n";
#  	}
#      }

#      for($i=0; $i<$nHits;$i++) {
#  	$DEBUG && print "DEBUG: $queryName $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitStrand[$i]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ($subjectHitStrand[$i]) ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]\n";  
#      }

#      # Then pop up a user interaction window for choosing hits for use in annotation.
#      my $qblast_win=$main->Toplevel;
#      $qblast_win->title("qblastn at NCBI results");
#      $qblast_win->geometry('+300+300');
#      $qblast_win->configure(-background=>'linen',-width=>'600');
    
#      my $qblast_main_frame=$qblast_win->Frame(-background=>$default_win_background)->pack(-fill => 'both', -expand=> 'yes');
#      my $qblast_list_frame=$qblast_main_frame->Frame(-background=>$default_win_background)->pack(-fill => 'both', -expand=> 'yes',-side=>'top');
#      my $qblast_hits_list=$qblast_list_frame->Listbox(-relief => 'sunken',-height => 25, -setgrid=>'true', -selectmode=> 'multiple')->pack(-expand=>'yes',-fill=>'both',-side=>'left');
#      for($i=0; $i<$nHits;$i++) {
#  	$hitentry="$queryName $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitStrand[$i]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ($subjectHitStrand[$i]) ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]";
#  	$DEBUG && print "DEBUG: $hitentry\n";
#  	$qblast_hits_list->insert('end',$hitentry);
#      }  

#      my $qblast_list_sby=$qblast_list_frame->Scrollbar(-command => ['yview', $qblast_hits_list])->pack(-side=>'right',-fill=>'y');
#      my $qblast_list_sbx=$qblast_main_frame->Scrollbar(-orient=>'horiz',-command => ['xview', $qblast_hits_list])->pack(-side=>'top',-fill=>'x');
#      $qblast_hits_list->configure(-yscrollcommand => ['set', $qblast_list_sby],-xscrollcommand => ['set', $qblast_list_sbx] );
#      # scrollbars... connect commands...
    
#      my $qblast_action_frame=$qblast_main_frame->Frame(-background=>$default_win_background)->pack(-fill => 'x', -expand=> 'yes',-side=>'top',-anchor=>'w');
#      my $qblast_annotate=$qblast_action_frame->Button(-command=>sub { # Q&D for finding our annotation id...
#  	#$annotation_id="ab_" . $queryName=~m/^[\w\d\.]+_(\d+)/;
#  	$annotation_nr=annotation_what_nr($annotation,$annotation_id);
#  	$note=$$annotation{note}->[$annotation_nr];
#  	foreach $selected ($qblast_hits_list->curselection) {
#  	    $note.="$queryName $queryHitBegin[$selected]-$queryHitEnd[$selected] ($queryHitStrand[$selected]) hit $subjectName[$selected] $subjectHitBegin[$selected]-$subjectHitEnd[$selected] ($subjectHitStrand[$selected]) ids=$identities[$selected]/$alignLength[$selected] score=$score[$selected] Expect=$expect[$selected]\n";
#  	}
#  	main::annotation_edit($main,$canvas,$annotation,$annotation_nr,$note);
#      },-text=>"Annotate query")->pack(-side=>'left');
#      my $qblast_manual=$qblast_action_frame->Button(-text=>"Add as blasthits",-command=>sub {
#  	foreach $selected ($qblast_hits_list->curselection) {
#  	    $annotation_nr=annotation_what_nr($annotation,$annotation_id);
#  	    add_blasthit($canvas,$annotation,$sheet,$seq,'blastn',$queryName,$$annotation{start}->[$annotation_nr]+$queryHitBegin[$selected],$$annotation{start}->[$annotation_nr]+$queryHitEnd[$selected],$queryHitStrand[$selected],$subjectName[$selected], $subjectHitBegin[$selected],$subjectHitEnd[$selected],$subjectHitStrand[$selected],$identities[$selected],$alignLength[$selected],$score[$selected],$expect[$selected],$alignment[$selected]);
#  	}
#  	  main::level_layout($annotation, 'blast'); 
#  	  main::redraw_annotations($canvas,$annotation,$sheet,$seq); # slightly ugly..

#      })->pack(-side=>'left');

#      my $qblast_viewalign=$qblast_action_frame->Button(-text=>"Display alignments",-command=>sub {
#  	# fancy view?
#  	$annotation_nr=annotation_what_nr($annotation,$annotation_id);

#  	print "Alignments (qblastn) for annotation $$annotation{id}->[$annotation_nr] on $$annotation{seqName}->[$annotation_nr].\n";

#  	foreach $selected ($qblast_hits_list->curselection) {
#  	    print "$queryName $queryHitBegin[$selected]-$queryHitEnd[$selected] ($queryHitStrand[$selected]) hit $subjectName[$selected] $subjectHitBegin[$selected]-$subjectHitEnd[$selected] ($subjectHitStrand[$selected]) ids=$identities[$selected]/$alignLength[$selected] score=$score[$selected] Expect=$expect[$selected]\n";
#  	    print "\n$alignment[$selected]\n";
#  	}
#      })->pack(-side=>'left');


#      my $qblast_cancel=$qblast_action_frame->Button(-text=>"Cancel",-command=>sub { $qblast_win->destroy; })->pack(-side=>'right');

#      # Choose classification viewer? (Hierarcical display with "DETAILED LIST" extracted from kinetoplastid gene nomenclature pages.)
#      # Suggest name accordning to Gen. nomenclature "whitepaper"?

#  }

#  sub parse_blastx_reply {
#    my $main=shift;
#    my $canvas=shift;
#    my $sheet=shift;
#    my $seq=shift;
#    my $annotation=shift;
#    my $annotation_id=shift;
#    my $ncbi_reply=shift;
  
#    # head
#    # ref, about query, 
#    # areamap
#    # Sequences producing significant alignments
#    # <PRE>
#    #  <a name = > </a><a href>id</a>name
#    # Length
#    # score, ids, strand  
#    # </PRE>
#    #<PRE> 
#    # Database:
  
#    # Ok, for a crude first attempt: remove HTML-tags. We lose information this way, 
#    # but on the other hand it is possible to use previously written WU-blast parsing code..

#    $$ncbi_reply=~s/\<a name \= \d+\>/\>/g; # Introduce > before hits, by substitution with the click-map name tag..
#    $$ncbi_reply=~s/\<.+?\>//g;
#    $DEBUG && print "DEBUG: $$ncbi_reply\n";

#    # IMPROVEMENT: Outdated parser version. Please introduce available code for long hit-names etc.

#    my $nHits=0;
#    my $inHit=0;
#    my $inDetail=0;
#    my $detailField=0;

#    my @subjectHitBegin;
#    my @subjectHitEnd;
#    my @subjectName;
#    my @queryHitBegin;
#    my @queryHitEnd;
#    my @queryHitFrame;
#    my $queryName;
#    my @identities;
#    my @alignLength;
#    my @score;
#    my @expect;

#    my $next_row_is_alignment = 0;
    
#    foreach $_ (split(/\n{1}/,$$ncbi_reply)) {    

#      # TAKE CARE WHEN EDITING ANY PATTERN - A FEW OF THEM OCCUR IN SEVERAL PLACES, SO CHANGE THEM ALL!
    
#      # Either we are looking at a detailed view of a hit, or we are scanning the hit-header.
#      # The rest of the output is kindly enough different.. =)
#      if($inHit) {
#        if($inDetail) { 
#  	  if($next_row_is_alignment) {
#  	      $alignment[$nHits-1] .= "$_\n";
#  	      $next_row_is_alignment = 0;
#  	  }
#  	if(/^Sbjct\:/) {
#  	  if($detailField==0) {
#  	    # Get both hitBegin and hitEnd on first encounter, then only the hitEnds..
#  	    ($subjectHitBegin[$nHits-1],$subjectHitEnd[$nHits-1])=/^Sbjct\:\s+(\d+)\s+[\w*-]+\s+(\d+)/;	    
#  	  } else { 
#  	    ($subjectHitEnd[$nHits-1])=/^Sbjct\:\s+\d+\s+[\w*-]+\s+(\d+)/;
#  	  }
#  	  $alignment[$nHits-1] .= "$_\n";
#  	} elsif(/^Query\:/) {
#  	  # A new block of hit alignment rows was found.
#  	  $detailField++;
#  	  ($queryHitEnd[$nHits-1])=/^Query\:\s+\d+\s+[\w*-]+\s+(\d+)/;
#  	  $alignment[$nHits-1] .= "\n$_\n";
#  	  $next_row_is_alignment = 1;
#  	} elsif(/Score\s{1}\=/) {
#  	  # Parse of alignment rows found a new hit on the same subject sequence.
#  	  $subjectName[$nHits]=$subjectName[$nHits-1];
#  	  $inDetail=0;
#  	  $nHits++;
#  	  ($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
#  #	  ($p[$nHits-1])=/P\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
#  	  ($expect[$nHits-1])=/Expect\s*\=\s+([0-9\.e\-\+]+)/;
#  	} elsif(/^\>/) {
#  	  # Parse of alignment rows found hit on new subject sequence.
#  	  $inDetail=0;
#  	  # Hits supposedly begin with a row "> FASTA_NAME"
#  	  ($subjectName[$nHits])=/^\>(.+)/;
#  	  $nHits++;
#  	  $alignment[$nHits-1]="";
#  	} 
#  	# in detail ends
#  	#Parse a hit header..     
#        } elsif(/Score\s{1}\=/) {
#  	($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
#  	# Syntax of the P value varies btw runs in blastn.. *sigh*
#  	($expect[$nHits-1])=/Expect\s+=\s+([0-9\.e\-\+]+)/;
#        } elsif (/Identities/) {
#  	($identities[$nHits-1],$alignLength[$nHits-1])=/Identities\s+\=\s+(\d+)\/(\d+)/;
#        } elsif (/Frame\s+\=/) {
#  	($queryHitFrame[$nHits-1])=/Frame\s+\=\s+([\d+-]+)/;
#        } elsif(/^Query\:/) {
#  	$inDetail=1;
#  	$detailField=0;
#  	# If this is a gapped alignment, the aligned sequences may contain dashes for gaps.. 
#  	($queryHitBegin[$nHits-1],$queryHitEnd[$nHits-1])=/^Query\:\s+(\d+)\s+[\w*-]+\s+(\d+)/;
#  	# Get both hitBegin and hitEnd on first encounter, later only the hitEnds..
#  	$alignment[$nHits-1] .= "$_\n";
#  	$next_row_is_alignment = 1;
#        }   
#      } elsif(/^\>/) {
#        # Hits supposedly begin with a row "> FASTA_NAME"
#        ($subjectName[$nHits])=/^\>(.+)/;
#        $nHits++;
#        $inHit=1;
#      } elsif(/^Query\=/) {
#        # Actually just evaluated once to get query name..
#        ($queryName)=/^Query\=\s*(.+)/; 
#  #      print "Query $queryName.\n";
#      }
#    }

#    for($i=0; $i<$nHits;$i++) {
#      $DEBUG && print "DEBUG: $queryName $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitFrame[$i]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]\n";  
#    }

#    # Then pop up a user interaction window for choosing hits for use in annotation.
#    my $qblast_win=$main->Toplevel;
#    $qblast_win->title("qblastx at NCBI results");
#    $qblast_win->geometry('+300+300');
#    $qblast_win->configure(-background=>'linen',-width=>'600');
  
#    my $qblast_main_frame=$qblast_win->Frame(-background=>$default_win_background)->pack(-fill => 'both', -expand=> 'yes');
#    my $qblast_list_frame=$qblast_main_frame->Frame(-background=>$default_win_background)->pack(-fill => 'both', -expand=> 'yes',-side=>'top');
#    my $qblast_hits_list=$qblast_list_frame->Listbox(-relief => 'sunken',-height => 25, -setgrid=>'true', -selectmode=> 'multiple')->pack(-expand=>'yes',-fill=>'both',-side=>'left');
#    for($i=0; $i<$nHits;$i++) {
#      $hitentry="$queryName $queryHitBegin[$i]-$queryHitEnd[$i] ($queryHitFrame[$i]) hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]";
#      $DEBUG && print "DEBUG: $hitentry\n";
#      $qblast_hits_list->insert('end',$hitentry);
#    }  

#    my $qblast_list_sby=$qblast_list_frame->Scrollbar(-command => ['yview', $qblast_hits_list])->pack(-side=>'right',-fill=>'y');
#    my $qblast_list_sbx=$qblast_main_frame->Scrollbar(-orient=>'horiz',-command => ['xview', $qblast_hits_list])->pack(-side=>'top',-fill=>'x');
#    $qblast_hits_list->configure(-yscrollcommand => ['set', $qblast_list_sby],-xscrollcommand => ['set', $qblast_list_sbx] );
#    # scrollbars... connect commands...
  
#    my $qblast_action_frame=$qblast_main_frame->Frame(-background=>$default_win_background)->pack(-fill => 'x', -expand=> 'yes',-side=>'top',-anchor=>'w');
#    my $qblast_annotate=$qblast_action_frame->Button(-command=>sub { # Q&D for finding our annotation id...
#  						     #$annotation_id="ab_" . $queryName=~m/^[\w\d\.]+_(\d+)/;
#  						     $annotation_nr=annotation_what_nr($annotation,$annotation_id);
#  						     $note=$$annotation{note}->[$annotation_nr];
#  						     foreach $selected ($qblast_hits_list->curselection) {
#  						       $note.="$queryName $queryHitBegin[$selected]-$queryHitEnd[$selected] ($queryHitFrame[$selected]) hit $subjectName[$selected] $subjectHitBegin[$selected]-$subjectHitEnd[$selected] ids=$identities[$selected]/$alignLength[$selected] score=$score[$selected]ex Expect=$expect[$selected]\n";
#  						     }
#  						     main::annotation_edit($main,$canvas,$annotation,$annotation_nr,$note);
#  						   },-text=>"Annotate query")->pack(-side=>'left');
#    my $qblast_manual=$qblast_action_frame->Button(-text=>"Add as blasthits",-command=>sub {
#        foreach $selected ($qblast_hits_list->curselection) {
#  	  $annotation_nr=annotation_what_nr($annotation,$annotation_id);
#  	  add_blasthit($canvas,$annotation,$sheet,$seq,'blastx',$queryName,$$annotation{start}->[$annotation_nr]+$queryHitBegin[$selected],$$annotation{start}->[$annotation_nr]+$queryHitEnd[$selected],$queryHitFrame[$selected],$subjectName[$selected], $subjectHitBegin[$selected],$subjectHitEnd[$selected],$identities[$selected],$alignLength[$selected],$score[$selected],$expect[$selected],$alignment[$selected]);
#        }
#  	  main::level_layout($annotation, 'blast'); 
#  	  main::redraw_annotations($canvas,$annotation,$sheet,$seq); # slightly ugly..
#    })->pack(-side=>'left');

#    my $qblast_viewalign=$qblast_action_frame->Button(-text=>"Display alignments",-command=>sub {
#        # fancy view?
#        $annotation_nr=annotation_what_nr($annotation,$annotation_id);

#        print "Alignments (qblastx) for annotation $$annotation{id}->[$annotation_nr] on $$annotation{seqName}->[$annotation_nr].\n";

#        foreach $selected ($qblast_hits_list->curselection) {
#  	  print "$queryName $queryHitBegin[$selected]-$queryHitEnd[$selected] ($queryHitFrame[$selected]) hit $subjectName[$selected] $subjectHitBegin[$selected]-$subjectHitEnd[$selected] ids=$identities[$selected]/$alignLength[$selected] score=$score[$selected] Expect=$expect[$selected]\n";
#  	  print "\n$alignment[$selected]\n";
#        }
#    })->pack(-side=>'left');

#    my $qblast_cancel=$qblast_action_frame->Button(-text=>"Cancel",-command=>sub { $qblast_win->destroy; })->pack(-side=>'right');

#    # Choose classification viewer? (Hierarcical display with "DETAILED LIST" extracted from kinetoplastid gene nomenclature pages.)
#    # Suggest name accordning to Gen. nomenclature "whitepaper"?

#  }

#  sub parse_blastp_reply {
#    my $main=shift;
#    my $canvas=shift;
#    my $sheet=shift;
#    my $seq=shift;
#    my $annotation=shift;
#    my $annotation_id=shift;
#    my $ncbi_reply=shift;
  
#    # head
#    # ref, about query, 
#    # areamap
#    # Sequences producing significant alignments
#    # <PRE>
#    #  <a name = > </a><a href>id</a>name
#    # Length
#    # score, ids, strand  
#    # </PRE>
#    #<PRE> 
#    # Database:
  
#    # Ok, for a crude first attempt: remove HTML-tags. We lose information this way, 
#    # but on the other hand it is possible to use previously written WU-blast parsing code..

#    $$ncbi_reply=~s/\<a name \= \d+\>/\>/g; # Introduce > before hits, by substitution with the click-map name tag..
#    $$ncbi_reply=~s/\<.+?\>//g;
#    $DEBUG && print "DEBUG: $$ncbi_reply\n";

#    # IMPROVEMENT: Outdated parser version. Please introduce available code for long hit-names etc.

#    my $nHits=0;
#    my $inHit=0;
#    my $inDetail=0;
#    my $detailField=0;

#    my @subjectHitBegin;
#    my @subjectHitEnd;
#    my @subjectName;
#    my @queryHitBegin;
#    my @queryHitEnd;
#    my $queryName;
#    my @identities;
#    my @alignLength;
#    my @score;
#    my @expect;

#    my $next_row_is_alignment = 0;
    
#    foreach $_ (split(/\n{1}/,$$ncbi_reply)) {    

#      # TAKE CARE WHEN EDITING ANY PATTERN - A FEW OF THEM OCCUR IN SEVERAL PLACES, SO CHANGE THEM ALL!
    
#      # Either we are looking at a detailed view of a hit, or we are scanning the hit-header.
#      # The rest of the output is kindly enough different.. =)
#      if($inHit) {
#        if($inDetail) { 
#  	  if($next_row_is_alignment) {
#  	      $alignment[$nHits-1] .= "$_\n";
#  	      $next_row_is_alignment = 0;
#  	  }
#  	if(/^Sbjct\:/) {
#  	  if($detailField==0) {
#  	    # Get both hitBegin and hitEnd on first encounter, then only the hitEnds..
#  	    ($subjectHitBegin[$nHits-1],$subjectHitEnd[$nHits-1])=/^Sbjct\:\s+(\d+)\s+[\w*-]+\s+(\d+)/;
#  	  } else { 
#  	    ($subjectHitEnd[$nHits-1])=/^Sbjct\:\s+\d+\s+[\w*-]+\s+(\d+)/;
#  	  }
#  	  $alignment[$nHits-1] .= "$_\n";
#  	} elsif(/^Query\:/) {
#  	  # A new block of hit alignment rows was found.
#  	  $detailField++;
#  	  ($queryHitEnd[$nHits-1])=/^Query\:\s+\d+\s+[\w*-]+\s+(\d+)/;
#  	  $alignment[$nHits-1] .= "\n$_\n";
#  	  $next_row_is_alignment = 1;
#  	} elsif(/Score\s{1}\=/) {
#  	  # Parse of alignment rows found a new hit on the same subject sequence.
#  	  $subjectName[$nHits]=$subjectName[$nHits-1];
#  	  $inDetail=0;
#  	  $nHits++;
#  	  ($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
#  #	  ($p[$nHits-1])=/P\(*\d*\)*\s+=\s+([0-9\.e\-\+]+)/;
#  	  ($expect[$nHits-1])=/Expect\s*\=\s+([0-9\.e\-\+]+)/;
#  	} elsif(/^\>/) {
#  	  # Parse of alignment rows found hit on new subject sequence.
#  	  $inDetail=0;
#  	  # Hits supposedly begin with a row "> FASTA_NAME"
#  	  ($subjectName[$nHits])=/^\>(.+)/;
#  	  $nHits++;
#  	  $alignment[$nHits-1]="";
#  	} 
#  	# in detail ends
#  	#Parse a hit header..     
#        } elsif(/Score\s{1}\=/) {
#  	($score[$nHits-1])=/Score\s+=\s+[\d\.]+\s*bits\s*\((\d+)\)/;
#  	# Syntax of the P value varies btw runs in blastn.. *sigh*
#  	($expect[$nHits-1])=/Expect\s+=\s+([0-9\.e\-\+]+)/;
#        } elsif (/Identities/) {
#  	($identities[$nHits-1],$alignLength[$nHits-1])=/Identities\s+\=\s+(\d+)\/(\d+)/;
#        } elsif(/^Query\:/) {
#  	$inDetail=1;
#  	$detailField=0;
#  	# If this is a gapped alignment, the aligned sequences may contain dashes for gaps.. 
#  	($queryHitBegin[$nHits-1],$queryHitEnd[$nHits-1])=/^Query\:\s+(\d+)\s+[\w*-]+\s+(\d+)/;
#  	# Get both hitBegin and hitEnd on first encounter, later only the hitEnds..
#  	$alignment[$nHits-1] .= "$_\n";
#  	$next_row_is_alignment = 1;
#        }   
#      } elsif(/^\>/) {
#        # Hits supposedly begin with a row "> FASTA_NAME"
#        ($subjectName[$nHits])=/^\>(.+)/;
#        $nHits++;
#        $inHit=1;
#      } elsif(/^Query\=/) {
#        # Actually just evaluated once to get query name..
#        ($queryName)=/^Query\=\s*(.+)/; 
#  #      print "Query $queryName.\n";
#      }
#    }

#    for($i=0; $i<$nHits;$i++) {
#      $DEBUG && print "DEBUG: $queryName $queryHitBegin[$i]-$queryHitEnd[$i] hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]\n";  
#    }

#    # Then pop up a user interaction window for choosing hits for use in annotation.
#    my $qblast_win=$main->Toplevel;
#    $qblast_win->title("qblastp at NCBI results");
#    $qblast_win->geometry('+300+300');
#    $qblast_win->configure(-background=>'linen',-width=>'600');
  
#    my $qblast_main_frame=$qblast_win->Frame(-background=>$default_win_background)->pack(-fill => 'both', -expand=> 'yes');
#    my $qblast_list_frame=$qblast_main_frame->Frame(-background=>$default_win_background)->pack(-fill => 'both', -expand=> 'yes',-side=>'top');
#    my $qblast_hits_list=$qblast_list_frame->Listbox(-relief => 'sunken',-height => 25, -setgrid=>'true', -selectmode=> 'multiple')->pack(-expand=>'yes',-fill=>'both',-side=>'left');
#    for($i=0; $i<$nHits;$i++) {
#      $hitentry="$queryName $queryHitBegin[$i]-$queryHitEnd[$i]  hit $subjectName[$i] $subjectHitBegin[$i]-$subjectHitEnd[$i] ids=$identities[$i]/$alignLength[$i] score=$score[$i] Expect=$expect[$i]";
#      $DEBUG && print "DEBUG: $hitentry\n";
#      $qblast_hits_list->insert('end',$hitentry);
#    }  

#    my $qblast_list_sby=$qblast_list_frame->Scrollbar(-command => ['yview', $qblast_hits_list])->pack(-side=>'right',-fill=>'y');
#    my $qblast_list_sbx=$qblast_main_frame->Scrollbar(-orient=>'horiz',-command => ['xview', $qblast_hits_list])->pack(-side=>'top',-fill=>'x');
#    $qblast_hits_list->configure(-yscrollcommand => ['set', $qblast_list_sby],-xscrollcommand => ['set', $qblast_list_sbx] );
#    # scrollbars... connect commands...
  
#    my $qblast_action_frame=$qblast_main_frame->Frame(-background=>$default_win_background)->pack(-fill => 'x', -expand=> 'yes',-side=>'top',-anchor=>'w');
#    my $qblast_annotate=$qblast_action_frame->Button(-command=>sub { # Q&D for finding our annotation id...
#  						     #$annotation_id="ab_" . $queryName=~m/^[\w\d\.]+_(\d+)/;
#  						     $annotation_nr=annotation_what_nr($annotation,$annotation_id);
#  						     $note=$$annotation{note}->[$annotation_nr];
#  						     foreach $selected ($qblast_hits_list->curselection) {
#  						       $note.="$queryName $queryHitBegin[$selected]-$queryHitEnd[$selected] hit $subjectName[$selected] $subjectHitBegin[$selected]-$subjectHitEnd[$selected] ids=$identities[$selected]/$alignLength[$selected] score=$score[$selected] Expect=$expect[$selected]\n";
#  						     }
#  						     main::annotation_edit($main,$canvas,$annotation,$annotation_nr,$note);
#  						   },-text=>"Annotate query")->pack(-side=>'left');

#    my $qblast_manual=$qblast_action_frame->Button(-text=>"Add as blasthits",-command=>sub {
#        foreach $selected ($qblast_hits_list->curselection) {
#  	  $annotation_nr=annotation_what_nr($annotation,$annotation_id);
#  	  add_blasthit($canvas,$annotation,$sheet,$seq,'blastp',$queryName,$$annotation{start}->[$annotation_nr]+$queryHitBegin[$selected],$$annotation{start}->[$annotation_nr]+$queryHitEnd[$selected],$$annotation{frame}->[$annotation_nr],$subjectName[$selected],$subjectHitBegin[$selected],$subjectHitEnd[$selected],$identities[$selected],$alignLength[$selected],$score[$selected],$expect[$selected],$alignment[$selected]);
#  	  $qblast_win->update();
#        }
#      main::level_layout($annotation, 'blast'); 
#      main::redraw_annotations($canvas,$annotation,$sheet,$seq); # slightly ugly..
      
#    })->pack(-side=>'left');

#    my $qblast_viewalign=$qblast_action_frame->Button(-text=>"Display alignments",-command=>sub {
#        # fancy view?
#        $annotation_nr=annotation_what_nr($annotation,$annotation_id);

#        print "Alignments (qblastp) for annotation $$annotation{id}->[$annotation_nr] on $$annotation{seqName}->[$annotation_nr].\n";

#        foreach $selected ($qblast_hits_list->curselection) {
#  	  print "$queryName $queryHitBegin[$selected]-$queryHitEnd[$selected] hit $subjectName[$selected] $subjectHitBegin[$selected]-$subjectHitEnd[$selected] ids=$identities[$selected]/$alignLength[$selected] score=$score[$selected] Expect=$expect[$selected]\n";
#  	  print "\n$alignment[$selected]\n";
#        }
#    })->pack(-side=>'left');

#    my $qblast_cancel=$qblast_action_frame->Button(-text=>"Cancel",-command=>sub { $qblast_win->destroy; })->pack(-side=>'right');

#    # Choose classification viewer? (Hierarcical display with "DETAILED LIST" extracted from kinetoplastid gene nomenclature pages.)
#    # Suggest name accordning to Gen. nomenclature "whitepaper"?

#  }
