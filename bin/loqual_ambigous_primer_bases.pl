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
my %ambipos;
my %rambipos;

while(my $pattern = <FUZZPATTERN>) {
    chomp $pattern;
    if($pattern =~ /\#/) {
	next;
    } else {
	my @call = ("fuzznuc", "-pattern", $pattern, "-mismatch", "0", "-complement", "-sequence", $sequence_file, "-outfile", "fuzznuc.search.$pattern.gff", "-rformat", "gff");
	(system(@call) == 0) or die "Failed system @call: $?\n";
	push @patterns, $pattern;
	while($pattern =~ m/[^ATGCatgc]{1}/g) {
	    my $ambi = pos($pattern)-1;
	    push @{$ambipos{$pattern}}, $ambi;
	    my $rpos = length($pattern) - 1 - $ambi; # 0-based pos on sense strand; length one too much: last pos is L-1, ne! 
	    push @{$rambipos{$pattern}}, $rpos;
	}
	$DEBUG && print "Pattern $pattern is ambigous in positions_0: ", join(" ", @{$ambipos{$pattern}}), ".\n";
    }
}
close FUZZPATTERN;

# read quality values
my %qual;
my %qual_name;

get_fasta_qual($qual_file);

# combine fuzzhit and ambipositions to loqual..

foreach my $pattern ( keys %ambipos ) {    
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

	    if($hit_strand eq '+') {
		foreach my $ambi (@{$ambipos{$pattern}}) {
		    my $edit_qual_at = $hit_start-1 + $ambi;
		    $DEBUG && print "edit $hit_sequence pos $edit_qual_at ($hit_start -1 + $ambi) old qual ", ${$qual{$hit_sequence}}[$edit_qual_at], "\n";
		    ${$qual{$hit_sequence}}[$edit_qual_at] = 0;
		}
	    } elsif($hit_strand eq '-') {
		foreach my $ambi (@{$rambipos{$pattern}}) {
		    my $edit_qual_at = $hit_start-1 + $ambi; # rambipos already set on + strand of pattern
		    $DEBUG && print "edit $hit_sequence pos $edit_qual_at ($hit_start -1 + $ambi) old qual ", ${$qual{$hit_sequence}}[$edit_qual_at], "\n";
		    ${$qual{$hit_sequence}}[$edit_qual_at] = 0;
		}
	    }
	}
    }
    close GFF;
}

# pretty write quality file

open QUALOUT, ">$qual_out_file";

for my $qual_entry (keys %qual_name) {
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
