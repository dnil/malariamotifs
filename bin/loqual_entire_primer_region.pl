#!/usr/bin/perl -w

my $WARNING = 1;
my $DEBUG = 1;

# cli
my $pattern_file = $ARGV[0];
my $sequence_file = $ARGV[1];
my $qual_file = $ARGV[2];
my $qual_out_file = $ARGV[3];

# get a plaintext file with IUPAC patterns
open FUZZPATTERN, $pattern_file;
my @patterns;

while(my $pattern = <FUZZPATTERN>) {
    chomp $pattern;
    if($pattern =~ /\#/) {
	next;
    } else {
	my @call = ("fuzznuc", "-pattern", $pattern, "-pmismatch", "2", "-complement", "-sequence", $sequence_file, "-outfile", "fuzznuc.search.$pattern.gff", "-rformat", "gff");
	(system(@call) == 0) or die "Failed system @call: $?\n";
	push @patterns, $pattern;
    }
}
close FUZZPATTERN;

# read quality values
my %qual;
my %qual_name;
my @qual_order;

get_fasta_qual($qual_file);

# combine fuzzhit and ambipositions to loqual..

foreach my $pattern ( @patterns ) {
    open GFF, "<fuzznuc.search.$pattern.gff";
    
    while(my $gffrow = <GFF>) {
	chomp $gffrow;
	if($gffrow =~ /\#/) {
	    next;
	} else {
	    # weve got a hit
	    my @r = split(/\s+/, $gffrow);
	    my $hit_sequence = $r[0];
	    my $hit_start = $r[3];
	    my $hit_end = $r[4];
	    my $hit_strand = $r[6];
	    
	    my $qualarr_len = @{$qual{$hit_sequence}};
	    my $hit_len = $hit_end - $hit_start +1;

	    my $pre_primer_limit = 10;

#	    $DEBUG && print STDERR "DEBUG: hit ($hit_sequence, $hit_start - $hit_end, $hit_strand) within primer limit $primer_limit.";
#	    $DEBUG && print STDERR "DEBUG: hit ($hit_sequence, $hit_start - $hit_end, $hit_strand) within primer limit ",$qualarr_len - $primer_limit + 1," since seq len is $qualarr_len and limit $primer_limit.\n";
	    if( $hit_end <= $hit_len + $pre_primer_limit ) {
#		$DEBUG && print STDERR "DEBUG: hit ($hit_sequence, $hit_start - $hit_end, $hit_strand) within primer limit $primer_limit.\n";
		for ( my $i = 0; $i < $hit_end ; $i++ ) {
		    ${$qual{$hit_sequence}}[$i] = 0;
		}
	    } elsif ( $hit_start >= $qualarr_len - ($hit_len + $pre_primer_limit) + 1) {
#		$DEBUG && print STDERR "DEBUG: hit ($hit_sequence, $hit_start - $hit_end, $hit_strand) within primer limit ",$qualarr_len - $primer_limit + 1," since seq len is $qualarr_len and limit $primer_limit.\n";
		for ( my $i = $hit_start-1; $i < $qualarr_len ; $i++ ) {
		    ${$qual{$hit_sequence}}[$i] = 0;
		}
	    }
	}
    }
    close GFF;
}

# pretty write quality file

open QUALOUT, ">$qual_out_file";

for my $qual_entry (@qual_order) { # (keys %qual_name) {
    print QUALOUT ">",$qual_name{$qual_entry}, "\n";

    my $qualarr_len = @{$qual{$qual_entry}};
    for (my $i = 0; $i < $qualarr_len-1; $i++) {
	print QUALOUT ${$qual{$qual_entry}}[$i];
	if ( ($i+1) % 50 == 0 ) {
	    print QUALOUT "\n";
	} else {
	    print QUALOUT " ";
	}
    }
    # print last value
    print QUALOUT ${$qual{$qual_entry}}[$qualarr_len-1]."\n";
}

# done
sub get_fasta_qual {
    my $nQual=0;

    my $fastaQualFileName=shift;
    open(QUALFILE, "<$fastaQualFileName") || die "Quality fasta input file $fastaQualFileName open failed\n";

    $WARNING && print "Reading quality values...";
    while(my $qr = <QUALFILE>) {
	chomp $qr; # remove newline
	$qr=~s/\s+$//; # remove trailing whitespace
	
	if($qr=~/^\>/ ) {
	    $nQual++;	    
	    ($significant_qual_name) = ($qr=~ /^>(\S+)/);
	    push @qual_order, $significant_qual_name;

	    if(exists($qual_name{$significant_qual_name})) {
		$WARNING && print STDERR "WARNING: duplicate qual file entry $significant_qual_name.\n";
	    }
	    ($qual_name{$significant_qual_name}) = ($qr=~/^\>(.+)/);
	    @{$qual{$significant_qual_name}}=();
	} else {
	    push @{$qual{$significant_qual_name}}, split(/\s+/,$qr);
	}
    }

    # Done processing quality file
    $WARNING && print "done.\n";
    $DEBUG && print "Read $nQual fasta qual file entries.\n";
    close(QUALFILE);
  
    return();
}
