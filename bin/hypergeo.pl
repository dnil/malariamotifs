#!/usr/bin/perl -w

my $M = 250; # se dominants in total aln (per each contig ONCE =)
my $N = 205; # co dito

my $p = $M / ($M + $N);
my $q = $N / ($M + $N);

while(my $row = <STDIN>) {
    chomp $row;

    if($row =~ m/SEfrac/) {	
	my @r = split /\s+/, $row;

	my $se = $r[1];
	my $co = $r[3];

	my $tries = $co + $se;
#	my $succ = $se;

	if( $tries == 0 ) {
	    print $row, "\tProb undefined.\n";
	    next;
	}
	
	my $hypgeo;
	if( $se != 0 ) {
	    $hypgeo = over($tries,$se)*  $p**$se * $q**$co ;
	} else {
	    $hypgeo = over($tries,$co)*  $p**$se * $q**$co ;
	}
	my $se_over = "-";
	if ($se/$tries >= $p) {
	    $se_over = "+";
	}
	
	print $row, "\tp ", sprintf("%.4f",$hypgeo)," ",$se_over,"\n";

    } else {

	print $row."\n";

    }
}


sub over {
    my $a = shift;
    my $b = shift;
    
    if ($b == 0) { 
	return(undef);
    }

    my $facb = fac($b);
    defined($facb) || return(undef);

    my $ret = pfac($a-$b+1,$a) / $facb ;
}

sub pfac {
    my $y = shift;
    my $x = shift;

    if($y > $x or $x < 1) {
	return(undef);
    }
    

    my $ret = 1;
    foreach my $i ($y..$x ) {
	$ret *= $i;
    }

    return($ret);
}


sub fac {
    my $x = shift;
    
    my $ret = 1;
    if($x > 1) {
	foreach my $i ( 2..$x ) {
	    $ret *= $i;
	}
    }
	    
    if($x == 0) {
	$ret = undef;
    }
    
    return($ret);
}
