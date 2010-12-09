#!/usr/bin/perl 

use Bio::SearchIO;

my $WARNING = 1;
my $DEBUG = 1;

#CLI
my $infile = $ARGV[0];

my $searchio = new Bio::SearchIO(-format => 'blast',
                                            '-file'   => $infile);


my %hittab;
my %flip;

while( my $result = $searchio->next_result ) {

    my $query_name = $result->query_name();
    $DEBUG && print "Parsing result for $query_name.\n";
    (defined($hittab{$query_name})) || ($hittab{$query_name} = {});
    $flip{$query_name} = -1; # ie no vote
    
    while( my $hit = $result->next_hit ) {
	my $query_length = $hit->query_length();
	my $hit_name = $hit->name();

	if($query_name eq $hit_name) { 
	    next; # ignore the hit to self (keep for debug!)
	}

	defined($hittab{$query_name}->{$hit_name}) || ($hittab->{$query_name}->{$hit_name} = 0);
	
	while( my $hsp = $hit->next_hsp ) {
	    if( $hsp->evalue < 1e-10 && $hsp->length / $query_length > 0.4 ) {

		($query_strand, $hit_strand) = $hsp->strand('list'); # rtfsc fix, but contrary to documentation.
		
		if ($query_strand eq $hit_strand) {
		    $DEBUG && print "Logging ++ ($query_strand $hit_strand) hit $query_name $hit_name.\n";
		    $hittab{$query_name}->{$hit_name} = 1;
		} elsif ($query_strand ne $hit_strand) {
		    $DEBUG && print "Logging +- ($query_strand $hit_strand) hit $query_name $hit_name.\n";
		    $hittab{$query_name}->{$hit_name} = -1;
		}
	    }
	}	
    }
}

my @recheck_stack;

foreach my $query_name ( keys %hittab ) {
    foreach my $hit_name ( keys %{$hittab{$query_name}} ) {

	if(  $hittab{$query_name}->{$hit_name} == 1 ) {
#	    $DEBUG && print "checking a ++ $query_name $hit_name.\n";
	    # check reciprocity? - shoudn't be neccessary

	    if( $flip{$query_name} == 1 ) {
		# ++ hit
		if($flip{$hit_name} == 1) {
		    # we have already flipped both hit and query, so a ++ hit is consistent with this
		} else { 
		    if($flip{$hit_name} == 0) {
			$DEBUG && $WARNING && print STDERR "WARNING: flipped query $query_name hits unflipped $hit_name ++, so we're going to flip $hit_name. But, an (other?) flipped & checked sequence has a +- hit to $query_name as well. Better check this one manually?\n";
			push @recheck_stack, $hit_name, $query_name;
		    }
		    # also if flip{hit_name} is unset (-1)
		    $flip{$hit_name} = 1;
		}

		
	    } else {
		# else, all is well: ++ and no previous  flip
		$DEBUG &&print "No flip on ++ $query_name $hit_name hit.\n";
		$flip{$query_name} = 0;
		$flip{$hit_name} = 0;
	    }
	} elsif( $hittab{$query_name}->{$hit_name} == -1 ) {
#	    $DEBUG && print "checking a +- $query_name $hit_name.\n";

	    if( $flip{$query_name} == 0 ) {
		if($flip{$hit_name} == 0) {
		    $DEBUG && $WARNING && print STDERR "WARNING: hit $query_name to $hit_name is +-, so we're going to flip $hit_name. But, an there is at least one other sequence with a hit driving $hit_name to remain unflipped as well. Better check this one manually?\n";
		    push @recheck_stack, $hit_name, $query_name;
		}
		# also if flip{hit_name} is unset (-1)
		$flip{$hit_name} = 1;
	    } elsif ($flip{$query_name} == 1)  { # else, all is well: +- and a previously registered flip. Lets register a non-flip for the pair for consistency checking.
		$DEBUG && print "No flip on +- $query_name $hit_name hit.\n";
		# note again: query_name already flipped.
		$flip{$hit_name} = 0;
	    } elsif ($flip{$query_name} == -1)  {
		# else, all is well: +- and a previously registered flip. Lets register a non-flip for the pair for consistency checking.
		$DEBUG && print "No flip on +- $query_name $hit_name hit.\n";
		# note again: query_name already flipped.
		$flip{$query_name = 0};
		$flip{$hit_name} = 0;
	    }
	}
    }
}

while (@recheck_stack) {
    $hit_name = shift @recheck_stack;
    $query_name = shift @recheck_stack;
    
    if(  $hittab{$query_name}->{$hit_name} == 1 ) {
	if( $flip{$query_name} == $flip{$hit_name} ) {
	    $DEBUG && $WARNING && print "recheck of $query_name to $hit_name ++ hit ok - equal flip status.\n";
	} else {
	    $WARNING && print "WARNING! Recheck of $query_name to $hit_name ++ hit: flip query is ",$flip{$query_name}, " but flip hit is ", $flip{$hit_name}, ". Please resolve manually - there should be more conflicting hits!\n";
	} 
    } elsif($hittab{$query_name}->{$hit_name} == -1) {
	if( $flip{$query_name} != $flip{$hit_name} ) {
	    $DEBUG && $WARNING && print "recheck of $query_name to $hit_name +- hit ok - different flip status.\n";
	} else {
	    $WARNING && print "WARNING! Recheck of $query_name to $hit_name +- hit: flip query is ",$flip{$query_name}, " but flip hit is ", $flip{$hit_name}, ". Please resolve manually - there should be more conflicting hits!\n";
	}
    }
}

# output list of sequences to flip
$DEBUG && print "Flip the following:\n";
foreach my $query_name ( keys %flip ) {
    if($flip{$query_name}) {
	print $query_name."\n";
    }
}
	
# a       -1  => flip d
# b       -1  => flip e
# c       0
# d=r(a)  -1  but flipped
# e=r(b)  -1

# let first hit decide global orientation
# 
