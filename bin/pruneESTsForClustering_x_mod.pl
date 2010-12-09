#!/usr/bin/perl -w
#
# Prune ESTs based on phred quality values
#
# Daniel Nilsson, 991201
# 991203 added code for multiple good stretches in one seq and hence rephrased the good stretch criteria.
# Removed that feature again for use in clustering and submission
# 
# Usage: 
# 
# pruneESTs seq_fasta_file qual_fasta_file pruned_ests_fasta_file [--prune-x [--write-qual pruned_qual_fasta_file]]
#

# Declare useful subs

sub getFastaSeq;
sub getFastaQual;
sub sum;

# Ajustable parameters

my $windowLen=15;
my $halfWindowLen=7;
my $goodLimit=16;
my $badLimit=17;
my $goodNumberOfBases=10;
my $badNumberOfBases=5;
my $shortestKeptGood=100;

my $nFASTA;
my @estSeq;
my @estName;
my $nQual;
my @estQual;

my @qualJudge;          # Array of judgement strings, i e [GBU].

# needless set of vectors.
my @goodies;            # Array of integers stating the number of good stretches in this sequence.

my @beginVectFree;
my @endVectFree;
my @beginGood;          # Array of good begins.
my @endGood;            # Array of happy endings.
my @badSeq;             # The bad ones. 

# Read input parameters

if((@ARGV<3) or ($ARGV[0] eq ("--help"|"-?"))) {
  die "Usage: pruneBacends.pl seq_fasta_file qual_fasta_file pruned_ests_fasta_file [--prune-x [--write-qual pruned_qual_fasta_file]]\n";
}

my $prunex=0;

if($ARGV[3] eq "--prune-x") {
  $prunex=1;
}

if($ARGV[4] eq "--write-qual") {
  $prunequal=1;
  $pruned_qual_fasta_file_name=$ARGV[5];
  ($pruned_qual_fasta_file_name eq $ARGV[1]) && die "Fatal error: Output file $pruned_qual_fasta_file_name same as input qual name.\n";
  open(PRUNEDQUAL, ">$pruned_qual_fasta_file_name") || die "Fatal error: Opening output file $pruned_qual_fasta_file_name failed.\n";
}

if($ARGV[2]) {
  open(PRUNEDEST, ">$ARGV[2]") || die "Output file $ARGV[2] open failed\n";
} else {
  print "No output file specified.\n";
  # or maby set PRUNEDEST=STDOUT? Hmm..
}; 

# First, get the sequences
#print "Reading sequences...\n";

@_=getFastaSeq($ARGV[0]);
$nFASTA=shift(@_);
@estSeq=@_[0..($nFASTA-1)];
@estName=@_[($nFASTA)..(2*$nFASTA-1)];

# $nr=@estSeq;
# print "Got $nr seqs ";
# $nr=@estName;
# print "and $nr names.\n"; 

@_=getFastaQual($ARGV[1]);
$nQual=shift(@_);
@estQual=@_;

# DEBUG
# $nr=@estQual;
# print "Got $nr quality value sequences.\n";

if($nQual != $nFASTA) {
  die "ERROR! Number of sequence entries does not equal number of quality entries.\n"
} 

print "Judging sequence quality...\n";

my $goods;
my $bads;

for($i=0;$i<$nQual;$i++) {

  # Fill the quality judgement string with U's,

  $qualJudge[$i]=$estSeq[$i]; 
  $qualJudge[$i]=~s/[atgcnxATGCNX]/U/g;
  
  $idiotic_comma = chop $estQual[$i];
  ($idiotic_comma eq ',') || die("Error in qual file handling. No end comma in ".$estName[$i].".\n");

  my @quals= split(/,+/,$estQual[$i]);
  
  # print "DEBUG: ($estName[$i]) @quals\n";

  # Decide ("judge") if a base is "G"ood or "B"ad based on two threshold values taken 
  # over a sliding window.
  #   
  #                  |----J----|                                          Sliding window
  #      ---------------------------------------------------------------  EST
  #
  # The first $halfWindowLen bases need special treatement, since a window centered on that base will 
  # be looking a $halfWindowLen outside the sequence.. (When window decides the central base judgement.)
  # My choice was to use the same window as before, but to judge the first base, 
  # depending on the overall sliding window sequence. This biases the quality of bases 
  # $halfWindowLen..$windowLen. It is btw already rather likely that they are Bad.
  #
  #  J---------|                                                          Sliding window
  # -------------------------------------------------------------         EST
  # 
  # In the same spirit, the final $halfWindowLen bases were judged by a window 
  # extending $windowLen backwards from the end..
  #                                                  |---------J          Sliding window
  # -------------------------------------------------------------         EST
  #

  # Initial part..

  for($j=0;$j<$halfWindowLen;$j++) {
    # Count good and bad bases in sliding window 
    $goods=0;
    $bads=0;       
    foreach(@quals[$j..($j+$windowLen)]) {      
      if ($_>$goodLimit) { 
	$goods++; 
      } elsif ($_<$badLimit) {
	$bads++;
      }
    }
    # and judge the leftmost base accordingly..
    $judgement="U";
    if ($goods>=$goodNumberOfBases) {
      $judgement="G"; 
    } elsif ($bads>=$badNumberOfBases) {
      $judgement="B";
    }
    # else they remain U.. 
    substr($qualJudge[$i],$j,1,$judgement);        
  } 
  
  # The central, normal, part of the EST..

  for($j=$halfWindowLen;$j<(@quals-$halfWindowLen);$j++) {   

    # Count good and bad bases in sliding window 
    $goods=0;
    $bads=0;  
    foreach (@quals[($j-$halfWindowLen)..($j+$halfWindowLen)]) {
      if ($_>$goodLimit) { 
	$goods++; 
      } elsif ($_<$badLimit) {
	$bads++;
      }
    }
    # and judge the window center accordingly..
    $judgement="U";
    if ($goods>=$goodNumberOfBases) {
      $judgement="G"; 
    } elsif ($bads>=$badNumberOfBases) {
      $judgement="B";
    }
    # else they remain U.. 

    substr($qualJudge[$i],$j,1,$judgement);    
  }

  # And then the same thing in the end; here the final base in the window will count, and 
  # bases @quals-$windowLen..@quals-$halfWindowLenth will bias the result..

  for($j=@quals-$halfWindowLen;$j<@quals;$j++) {

    # Count good and bad bases in sliding window 
    $goods=0;
    $bads=0;  
    foreach(@quals[scalar($j-$windowLen)..$j]) {
      if ($_>$goodLimit) { 
	$goods++; 
      } elsif ($_<$badLimit) {
	$bads++;
      }
    }
    # and judge the right-most base accordingly..
    $judgement="U";
    if ($goods>=$goodNumberOfBases) {
      $judgement="G"; 
    } elsif ($bads>=$badNumberOfBases) {
      $judgement="B";
    }
    # else they remain U.. 
    substr($qualJudge[$i],$j,1,$judgement);
  } 

  # Trivial pruning: retain only the interval [first base judged good .. first base judged bad]
  # NB. THAT FAILS MISERABLY ON MULTIPLE GOOD SUBSEQUENCES.. SPLIT them if longer than $shortestKeptGood, 
  # otherwise dump'em.
  #
  # Ran into trouble, since it _seems_ any ^ or $ in a pattern resets the search engine, so a scalar m// 
  # won't do the job. Implementing a little search machine of my one using "index" and some support 
  # wariables.

  my $goOn=1;             # State machine state variable. 
  my $currentSearchPos=0; # Dito.
 
  $goodies[$i]=0;         # No goodies for free.
  
  #while($goOn) {  
  #
 
  # CLUSTER: For the purpose of clustering, get the longest possible "G"-stretch,
  # regardless of whether it contains bad bases or not. 

  # print "DEBUG: $qualJudge[$i]\n";
  
  $qualJudge[$i]=~m/B*(G+\w*G+)B*/;
  if($1) {
    $bestHit=$1;
  } else {
    $bestHit="";
  }

  # Still, dump sequence if the longest good sequence was too short.
  if(length($bestHit)>$shortestKeptGood) {
    $beginGood[$i]=index($qualJudge[$i],$bestHit);
    $endGood[$i]=$beginGood[$i]+length($bestHit) - 1; # 4b0, l 4 -> 8b0  # beginGood and endGood are inclusive coords
    $goodies[$i]=1;
  } else {
    push(@badSeq,$i);
#    print "$estName[$i] was found to have a too short best hi-qual stretch (".length($bestHit)." nt).\n"; # DEBUG?
    $goodies[$i]=0;
    next;
  }

  #    $beginGood[$i][$goodies[$i]]=index($qualJudge[$i],"G",$currentSearchPos);
  #    if($beginGood[$i][$goodies[$i]]==-1) {
  #      if ($goodies[$i]==0) {
  #	# If no good bases were found, this is a bad sequence. (No shit, Sherlock!)      
  #	push(@badSeq,$i);
  #	# So, deal with the next sequence.
  #	$goOn=0;
  #       } else {
  # 	# Ok, nothing more to see in here, folks. 
  # 	$goOn=0;
  #       }
  #     } else {
  #       # Looks like we may have found us some good sequence. Update search position.
  #       $currentSearchPos=$beginGood[$i][$goodies[$i]];
  #       # Well, now, where does this goodie end?
  #       $endGood[$i][$goodies[$i]]=index($qualJudge[$i],"B",$currentSearchPos);
  #       if ($endGood[$i][$goodies[$i]]==-1) {
  # 	# Thats as long as it gets, dude! (No more bad ones found!)
  # 	$endGood[$i][$goodies[$i]]=length($qualJudge[$i])-1;
  # 	# It's time to move on.
  # 	$goOn=0;
  #       } else {
  # 	# Ok, there it ends. We'll want to check the rest of the sequence to 
  # 	# see if there are more good stretches here, so update the search-position.
  # 	$currentSearchPos=$endGood[$i][$goodies[$i]]+1;
  #       }
  #       if(($endGood[$i][$goodies[$i]]-$beginGood[$i][$goodies[$i]])>$shortestKeptGood) {
  # 	# Save a later check - if the match was shorter than our shortestKeptGood, dont bother 
  # 	# counting it.
  
  # 	# It would seem that we have found a Good match! 
  # 	$goodies[$i]++;
  #       }
  #    }
  #  }

  # If good bases were found, set name of sequence to reflect what bases where judged 
  # good, so they can be retained.
  # For other purpouses than tests, we might now want to consider pruning it here..

  # no use pruning a discarded sequence further..
  if($goodies[$i]) {
      if($prunex) {
	  $beginVectFree[$i] = $beginGood[$i];
	  $endVectFree[$i] = $endGood[$i];
	  
	  $currentQualPrunedSeq = substr($estSeq[$i], $beginGood[$i], $endGood[$i] - $beginGood[$i] +1);
	  
	  ($startx)=($currentQualPrunedSeq=~m/^(X+)/);
	  ($endx)=($currentQualPrunedSeq=~m/(X+)$/);  

	  # drop sequence on multiple X-clusters
	  (@xclusters) = ($currentQualPrunedSeq=~m/[^X]+(X+)(?=[^X]+)/g);
	  if (@xclusters >1) {
	      print "Sequence $estName[$i] was found to have multiple vector clusters during vectorpruning and hence discarded.\n"; # DEBUG?
	      $goodies[$i]=0; # Throw it away! 
	      push(@badSeq,$i);
	      next;
	  }
	  
#	  ($middlex)=($currentQualPrunedSeq=~m/[^X]+(X+)[^X]+/);
	  $middlex = $xclusters[0]; 

	  if($startx) {
	      # print "DEBUG: initial X stretch: $startx\n";
	      $beginVectFree[$i] = $beginGood[$i] + length($startx); # length is 1 or more, but we want first _free_ base
	  }

	  if($endx) {
	      # print "DEBUG: final X stretch: $endx\n";
	      $endVectFree[$i] = $beginGood[$i] + length($currentQualPrunedSeq) -length($endx) -1;   # Last free base: 
	  }

	  if($middlex) {
	      # count total number of x:es? or never mind, since thats going to involve reading through all seq anyways, and do that after actual pruning?

#	      my @current_qp_seq = split(/ +/, $currentQualPrunedSeq);
	      if( $endGood[$i] - $beginGood[$i] +1 - length($middlex) < $shortestKeptGood) {
		  print "Sequence $estName[$i] was found too short during vectorpruning and discarded.\n"; # DEBUG?
		  $goodies[$i]=0; # Throw it away!
		  push(@badSeq,$i);
	      }

	      ($left, $right) = ( $currentQualPrunedSeq=~m/([^X]+)X+([^X]+)/ );

	      if (length($left) > length($right) ) {
		  if (length($left) > $shortestKeptGood) {
		      # print "DEBUG: initial X stretch: $startx\n";
                      # keeping left; moving end
		      $endVectFree[$i] = $beginGood[$i] + length($currentQualPrunedSeq)-1 -length($middlex) -length($right);  # length is 1 or more, but we want first _free_ base
		  } else {
		      print "Sequence $estName[$i] was found too short during vectorpruning and discarded (right side > left, but shorter than $shortestKeptGood).\n"; # DEBUG?
		      $goodies[$i]=0; # Throw it away!
		      push(@badSeq,$i);
		  }
	      } else {
		  if (length($right) > $shortestKeptGood) {
		      # keeping right; move start
		      $beginVectFree[$i] = $beginGood[$i]+length($left)+length($middlex);
		  } else {
		      print "Sequence $estName[$i] was found too short during vectorpruning and discarded (right side > left, but shorter than $shortestKeptGood).\n"; # DEBUG?
		      $goodies[$i]=0; # Throw it away! 
		      push(@badSeq,$i);
		  }
	      }
	  }

	  # print "DEBUG: $estName[$i]\n$currentQualPrunedSeq\n";
	  # Special case where all of the sequence is vector
	  if( $currentQualPrunedSeq=~m/^X+$/ ) {
	      print "$estName[$i] discarded as pure vector.\n"; # DEBUG?
	      $goodies[$i]=0; # Throw it away! 
	      push(@badSeq,$i);
	  }
      }
  }

  # Only zero or one goodie!
  if($goodies[$i]) {
      # Are we pruning away masked sequence as well?
      if($prunex) {
	  # Now, update beginGood and endGood if vect. free happen to more restrictive.
	  if($beginVectFree[$i]>$beginGood[$i]) {
	      $beginGood[$i]=$beginVectFree[$i];
	  }
	  if($endVectFree[$i]<$endGood[$i]) {
	      $endGood[$i]=$endVectFree[$i];
	  }
	  
	  if($endGood[$i]-$beginGood[$i]+1 < $shortestKeptGood) {
	      print "Sequence $estName[$i] was found too short after vectorpruning. Discarded.\n"; # DEBUG?
	      $goodies[$i]=0; # Throw it away! 
	      push(@badSeq,$i);
	  }	  
      }
  }      

  if($goodies[$i]) {
      print PRUNEDEST ">$estName[$i] GOOD: ",$beginGood[$i]+1,"-",$endGood[$i]+1;
#    print PRUNEDEST " VFree: ",$beginVectFree+1,"-",$endVectFree+1;
      print PRUNEDEST "\n";
      print PRUNEDEST substr($estSeq[$i],$beginGood[$i],$endGood[$i]-$beginGood[$i]+1),"\n";
      
      if($prunequal) {
	  print PRUNEDQUAL ">$estName[$i] GOOD: ",$beginGood[$i]+1,"-",$endGood[$i]+1;
	  print PRUNEDQUAL "\n";
#	  my @writequals=split(/,+/,$estQual[$i]);
	  my $qual_string = join(' ',@quals[$beginGood[$i]..$endGood[$i]]);
	  print ">$estName[$i] GOOD: ",$beginGood[$i]+1,"-",$endGood[$i]+1," VFree: ",$beginVectFree[$i]+1,"-",$endVectFree[$i]+1," length of current qualarr ", scalar(@quals)," seq ",length($estSeq[$i]),"\n"; # DEBUG
	  print PRUNEDQUAL $qual_string,"\n";
      }

  }
  
  # Check for left over U's.  
  if($qualJudge[$i]=~/U/) {
    print "WARNING: in pruneESTs.pl; sequence $i ($estName[$i]) had bases that fell between the good and bad tresholds.\n";
  }

  # Trial: print the judgements..
  #   print "The qualities of the sequence $estName[$i] were judged:\n";
  #   print $qualJudge[$i],"\n";
  #   if($goodies[$i]>0) {
      #     print "and pruneEST retained the following bases: ";
  #     print "$beginGood[$i]-$endGood[$i]";
  #     print ".\n";
  #   }
}

print "The following sequences were marked as trash:\n";
foreach(@badSeq) {
  print "Number $_: $estName[$_]\n"; 
}

################################# MAIN END ########################################

sub getFastaSeq {
  my $fastaFileName=shift;
  open(FASTAFILE,$fastaFileName ) || die "Sequence fasta input file $fastaFileName open failed.\n";

  my @fasta;
  my @name;

  # First, get the sequences
  print "Reading sequences...\n";
  my $nFASTA=0;
  while(<FASTAFILE>) {
    chomp;
    if(/^\>/) {
      # On fasta description line
      $nFASTA++;
      ($name[$nFASTA-1])=/^\>(.+)/;    
      $fasta[$nFASTA-1]="";
    } else {
      # Well, either the input is broken, or this is sequence data. Let us assume a sentient user.. :)
      # Get all genomic sequence chars into that $fasta string..
      s/[^atgcnxATGCNX]//g;
      s/\s+//g; # remove any pesky trailing whitespaces
      # UPPERCASE the "x"es
      tr/x/X/;

      $fasta[$nFASTA-1] .= $_;
    }
  }

  # Done processing fasta sequence file
  close(FASTAFILE);

  return ($nFASTA,@fasta,@name); 
}

sub getFastaQual {
  my $nQual=0;

  my $fastaQualFileName=shift;
  open(QUALFILE, "<$fastaQualFileName") || die "Quality fasta input file $fastaQualFileName open failed\n";

  print "Reading quality values...\n";
  while(my $qr = <QUALFILE>) {
    chomp $qr; # remove newline
    $qr=~s/\s+$//; # remove trailing whitespace

    if($qr =~ /^\>/) {
      # On description line.. Name field should be pretty much equal to the fasta file, so ignore it.	
      $nQual++;
      $qual[$nQual-1]="";
    } else {
	# beware of the blanklines
	if($qr ne "") {
	    $qual[$nQual-1].=join(',',split(/\s+/,$qr));
	    $qual[$nQual-1].=",";
	}
    }  
  }
    
  # Done processing quality file
  close(QUALFILE);
  
  return($nQual,@qual);
}

sub sum {
  $sum=0;
  # The pop-version looked nice, but choked on qual==0...
  for($n=0;$n<@_;$n++) {
    $sum+=$_[$n];
  }
  return $sum;
}
