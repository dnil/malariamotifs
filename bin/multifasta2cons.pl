#!/usr/bin/perl -w

use Bio::AlignIO;
use Bio::SeqIO;

my $inputfilename = $ARGV[0];
my $outputfilename = $ARGV[1];

my $in  = Bio::AlignIO->new(-file => $inputfilename ,
			 '-format' => 'fasta');

my $out = Bio::SeqIO->new(-file => ">$outputfilename" ,
			 '-format' => 'Fasta');


my ($outbasename) = ( $outputfilename =~ m{([^/]+)$} );

while ( my $aln = $in->next_aln() ) {

#    my $cons_str = consensus_string($aln, 26); # local thickness, so no limit on threshold
    my $cons_str = $aln->consensus_string(16);
    
#     if($aln->no_sequences > 100) {
# 	$cons_str = $aln->consensus_string(5); # trims away ends of thick alignments, will insert incorrect gaps regularly in the perhaps 100 / t thick sequences (thickness 5 at 20, 10 at 10..)
#     } elsif($aln->no_sequences > 50) {
# 	$cons_str = $aln->consensus_string(6); # trims away ends of thick alignments, will insert incorrect gaps regularly in the perhaps 100 / t thick sequences (thickness 5 at 20, 10 at 10..)
#     } elsif($aln->no_sequences > 10) {
# 	$cons_str = $aln->consensus_string(15); # trims away ends of thick alignments, will insert incorrect gaps regularly in the perhaps 100 / t thick sequences (thickness 5 at 20, 10 at 10..)
#     }  elsif($aln->no_sequences > 6) {
# 	$cons_str = $aln->consensus_string(16); # trims away ends of thick alignments, will insert incorrect gaps regularly in the perhaps 100 / t thick sequences (thickness 5 at 20, 10 at 10..)
#     } else {
# 	$cons_str = $aln->consensus_string(50); # trims away ends of thick alignments, will insert incorrect gaps regularly in the perhaps 100 / t thick sequences (thickness 5 at 20, 10 at 10..)

    $cons_str =~ s/\?//g; # remove ambigous positions (majority base had less than threshold votes)
    my $cons = Bio::Seq->new( '-seq' => $cons_str, '-id' => $outbasename );

    $out->write_seq( $cons );
}
		      
sub consensus_string {
    my $aln = shift;
    my $threshold = shift;

    my $out = "";
    my $len = $aln->length - 1;

#    foreach ( 0 .. $len ) {
	
#	$non_gap_count = ;
#	$smoothed_non_gap_count = ;
	
#    }

    foreach ( 0 .. $len ) {
        $out .= _consensus_aa($aln,$_,$threshold);
    }
    return $out;
}

sub _consensus_aa {
    my $aln = shift;
    my $point = shift;
    my $threshold_percent = shift || -1 ;
    my ($seq,%hash,$count,$letter,$key);
    my $gapchar = $aln->gap_char;

    my $gapcount;

    foreach $seq ( $aln->each_seq() ) {
        $letter = substr($seq->seq,$point,1);
	
        if ($letter eq $gapchar || $letter =~ /\./)  {
	    $gapcount++;
	    next;
	}
        # print "Looking at $letter\n";
        $hash{$letter}++;
    }

    my $number_of_sequences = $aln->no_sequences();
    my $number_of_voters = $number_of_sequences - $gapcount;
    my $min_voters = $number_of_voters * $threshold_percent / 100. ;
    $count = -1;
    $letter = '?';

    my $min_trusted_thickness = 6;
    my $trusted_thickness = ($treshold_percent * $number_of_sequences > $min_trusted_thickness) ? $treshold_percent * $number_of_sequences : $min_trusted_thickness ; # max (6, N% of total thickness) gives a sturdy, nonrandom coverage for the bases called by (tresholed) majority vote
    if ( $number_of_sequences < $min_trusted_thickness ) { 
	# very thin contig. Not much to do here without quality values but accept the majority decision
    } elsif ( $number_of_voters < $trusted_thickness ) {
	$letter = '?';
	return $letter;
    }


#    if( $number_of_voters > $trusted_thickness ) { 
#	
#    }

    foreach $key ( sort keys %hash ) {
        # print "Now at $key $hash{$key}\n";
        if( $hash{$key} > $count && $hash{$key} >= $threshold_percent ) {
            $letter = $key;
            $count = $hash{$key};
        }
    }

    return $letter;
}

