#!/usr/bin/perl -w

# /count_fractions_in_ace_file.pl < se31reads.fasta.ace | sort -k2,2nr -k5,5nr

my $isolate_sort = 0;
my $DEBUG = 0;
my $WARNING = 1;

while (my $arg = shift @ARGV) {
    if($arg eq "--debug") {
	$DEBUG = 1;
    } elsif($arg eq "--isolatefrac") { # list type of isolates in contig (se/uks/co/cp.. =)
	$isolate_sort = 1;
    } elsif ($arg eq "--no-warn") {	
	$WARNING = 0;
    } else {
	print STDERR "WARNING: Unrecognised command line parameter $arg. Continuing anyway.\n";
    }
}

my $current_contig ="";
my %contig = ();
my %contig_primer_count = ();
my $contig_members = -1;
my %clones_per_isolate = ();
my %clones_per_isolate_and_primer = ();

while(<STDIN>) {
    if (/^CO\s+(\S+)\s+\d+\s+(\d+)/) {
	$current_contig = $1;
	$contig_members = $2;

	$DEBUG && print "Contig $current_contig has $contig_members members.\n";
    }

    if(/^AF\s+(\S+)\s+/) {
        # count nr of forw/rev pairs, single rev or single forw

	my $read_name = lc $1;
	my ($direction,$isolate_type,$isolate_nr,$primer,$seqnr);

	if($read_name =~ /^([fr])([a-z]+)(\d+)([a-z]+)(\d+)/) {
	    ($direction,$isolate_type,$isolate_nr,$primer,$seqnr) = ($1,$2,$3,$4,$5);
	} else { 
	    $read_name =~ /^(\w+)(\d+)([a-z]+)(\d+)/;
	    $direction = "r";
	    ($isolate_type,$isolate_nr,$primer,$seqnr) = ($1,$2,$3,$4);
	    $WARNING && print STDERR "Found read $read_name without direction : setting reverse dir.\n";
	}

	$clone = $isolate_type.$isolate_nr.$primer.$seqnr;
	$DEBUG && print "$read_name\t$clone\n";
	
	if( !exists($contig{$clone})) {
	    $contig{$clone} = $current_contig;
	    
	    if( exists($contig_primer_count{$current_contig}) && exists(${$contig_primer_count{$current_contig}}->{$primer} ) ) {
		${$contig_primer_count{$current_contig}}->{$primer}++;
	    } else {
# 		$contig_primer_count{$current_contig} = \();
		${$contig_primer_count{$current_contig}}->{$primer} = 1;
	    }
		
	    if($isolate_sort) {
		if(exists(${$contig_isolate_type{$current_contig}}->{$isolate_type})) {
		    ${$contig_isolate_type{$current_contig}}->{$isolate_type}++;
		} else {
		    ${$contig_isolate_type{$current_contig}}->{$isolate_type} = 1;
		}
	    }

	    my $isolate = $isolate_type . $isolate_nr;
	    
	    if( exists(${${$contig_isolate_count{$current_contig}}->{$isolate}}->{$primer}) ) {
		${${$contig_isolate_count{$current_contig}}->{$isolate}}->{$primer}++;
	    } else {
		${${$contig_isolate_count{$current_contig}}->{$isolate}}->{$primer} = 1;
	    }

	    if( exists($clones_per_isolate{$isolate}) ) {
		$clones_per_isolate{$isolate}++;
	    } else {
		$clones_per_isolate{$isolate} = 1;
	    }

	    if( exists(${$clones_per_isolate_and_primer{$isolate}}->{$primer}) ) {
		${$clones_per_isolate_and_primer{$isolate}}->{$primer}++;
	    } else {
		${$clones_per_isolate_and_primer{$isolate}}->{$primer} = 1;
	    }

	} else {
	    $DEBUG && print "Ignoring $read_name since clone $clone already exists.\n";
	    # ignoring second found f/r for this clone
#	    ${$contig_primer_count{$current_contig}}->{$primer}++;
	}
    }        
}

#foreach my $isolate ( keys %clones_per_isolate ) {
   
#   $fraction_of
#   $isolate ) {


#}


foreach my $cont (keys %contig_primer_count) {
    
    my $number_of_clones_in_contig = 0;

    foreach my $primer (keys %${$contig_primer_count{$cont}}) {
	$DEBUG && print "Counting clones for primer $primer\n";
	$number_of_clones_in_contig += ${$contig_primer_count{$cont}}->{$primer};
    }
    
    if($isolate_sort) {
	foreach my $isolate_type ( keys %${$contig_isolate_type{$cont}}) { 	    
	    my $fraction_of_isolate_type_of_clones_in_contig = ${$contig_isolate_type{$cont}}->{$isolate_type} / $number_of_clones_in_contig;
	    print "$cont\t$isolate_type\t",sprintf("%.2f",$fraction_of_isolate_type_of_clones_in_contig),"\n";
	}
    } else {

	my %number_of_clones_in_isolate_in_contig;

	foreach my $isolate ( keys %${ $contig_isolate_count{$cont}} ) {
	    foreach my $primer ( keys %${${$contig_isolate_count{$cont}}->{$isolate}} ) {
		$DEBUG && print "Counting clones for isolate $isolate primer $primer\n";
		
		$primer_count = ${${$contig_isolate_count{$cont}}->{$isolate}}->{$primer};
		$number_of_clones_in_isolate_in_contig{$isolate} += $primer_count;
		
	    }
	}
    
#    my $fraction_of_primer_clones_in_isolate_in_contig;
	foreach my $isolate ( keys %${ $contig_isolate_count{$cont}} ) {
	    foreach my $primer ( keys %${${$contig_isolate_count{$cont}}->{$isolate}} ) {

#	    $fraction_of_clones_in_isolate_in_contig_of_primer = ${${${$contig_isolate_count{$cont}}->{$isolate}}}->{$primer}/$number_of_clones_in_isolate_in_contig{$isolate};	    
	    
		$fraction_of_clones_in_isolate_in_contig_out_of_total_isolate_count = $number_of_clones_in_isolate_in_contig{$isolate} 
		/ $clones_per_isolate{$isolate};	    
		
		my $fraction_of_clones_of_primer_in_contig_of_total_nr_of_primers_in_isolate = 
		    ${${$contig_isolate_count{$cont}}->{$isolate}}->{$primer} / ${$clones_per_isolate_and_primer{$isolate}}->{$primer};

		$primer_count = ${${$contig_isolate_count{$cont}}->{$isolate}}->{$primer};

		print $cont."\t".$number_of_clones_in_contig."\t".$isolate."\t".$number_of_clones_in_isolate_in_contig{$isolate}."\t".$clones_per_isolate{$isolate}."\t".$primer."\t".$primer_count."\t"
		    .sprintf("%.2f",$fraction_of_clones_in_isolate_in_contig_out_of_total_isolate_count)."\t"
		    .sprintf("%.2f",$fraction_of_clones_of_primer_in_contig_of_total_nr_of_primers_in_isolate)."\n";
	    }
	}
    }
#    my $fraction_of_primer_clones_in_contig;
#    foreach my $primer (keys %${$contig_primer_count{$cont}}) {
#	my $primer_count = ${$contig_primer_count{$cont}}->{$primer}; 
#	$fraction_of_primer_clones_in_contig = $primer_count / $number_of_clones_in_contig;
#    }

}

#foreach $clone (keys %contig) {  
#}
