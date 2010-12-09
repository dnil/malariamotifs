#!/usr/bin/perl -w
#
# ace2raw.pl 
#
# March 2005, Daniel Nilsson
#
# Usage: ace2raw.pl --infile <filename> [--outfile <filename>][--nowarn][]

use Bio::Assembly::IO;
use Bio::Assembly::Scaffold;
use Bio::Assembly::ScaffoldI;
use Bio::Assembly::Contig;

sub usage {
    print STDERR "Usage: ace2raw.pl --infile <filename> [--nowarn]\n";
    exit -1;
}

# parse command line arguments
my $DEBUG = 1;
my $WARNING = 1;

my $inFileName = undef;

while (my $arg = shift @ARGV) {
    if($arg eq "--debug") {
	$DEBUG = 1;
    } elsif ($arg eq "--infile") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No infile name given.");
	} else {
	    $inFileName = $arg;
	}
    } elsif ($arg eq "--nowarn") {
	$WARNING = 0;
    } else {
	usage();
    }
}

my $infile = Bio::Assembly::IO->new(-file=>"<$inFileName", -format=>"ace");
# my $outfile = Bio::Assembly::IO->new(-file=>">$outFileName", -format=>"raw");

$DEBUG && print STDERR "next assembly!\n";
my $assembly = $infile->next_assembly();

print STDERR "id:\t\t", $assembly->id, "\ncontigs:\t\t", scalar( @{[ $assembly->get_contig_ids() ]} ), "\nseqs in contigs:\t", $assembly->get_nof_sequences_in_contigs, "\n\n";
print STDERR "And the contigs are :", join(" ",$assembly->get_contig_ids()), "\n";

# print raw..

foreach my $contig ( $assembly->all_contigs() ) 
{    
    my $contig_id = $contig->id;
    open CONTIGOUT, ">$contig_id.raw";
    open CONTIGIDSOUT, ">$contig_id.ids";
    
    foreach my $seq ( $contig->each_seq() ) {
	# initgap!
	my $coord = $contig->get_seq_coord($seq);
	my $padded_start = $coord->start;

	my $initgap = ' ' x ($padded_start-1);
	print CONTIGOUT $initgap, $seq->seq(), "\n";
	print CONTIGIDSOUT $seq->id(), "\n";
    }
}
