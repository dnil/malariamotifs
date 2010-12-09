#!/usr/bin/perl -w
#
# fastare.pl 
#
# (c) August 2002 Daniel Nilsson
#
# apply iupac-aware perl regexp to nucleotide multifastafile
#
# Usage: fastare.pl [--infile <filename>][--outfile <filename>][--regexp <perl re>][--iupac][--case-insensitive][--no-revcomp]

# parse command line arguments
my $DEBUG = 0;

my $infh = *STDIN;
my $outfh = *STDOUT;

my $seqFileName = undef;
my $outFileName = undef;
my $re = undef;
my $iupac = 0;
my $caseins = 0;
my $revcomp = 1;

while (my $arg = shift @ARGV) {
    if($arg eq "--debug") {
	$DEBUG = 1;
    } elsif ($arg eq "--infile") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No infile name given.");
	} else {
	    $seqFileName = $arg;
	}
    } elsif ($arg eq "--outfile") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No outfile name given.");
	} else {
	    $seqFileName = $arg;
	}
    } elsif ($arg eq "--regexp") {
	$arg = shift @ARGV;
	if($arg eq "") {
	    usage("Error! No regular expression given.");
	} else {
	    $re = $arg;
	}
    } elsif ($arg eq "--case-insensitive") {
	$caseins = 1;
    } else {
	usage();
    }
}

!defined($re) && usage("Error! No regular expression given.");
$caseins && ($re = "(?i)$re");

# open infile or STDIN 
if(!defined($seqFileName)) {
    $DEBUG && debug("Reading input from STDIN."); # default behaviour
} else {
    local *INFILE;
    -e $seqFileName or fatal("No such infile $seqFileName.");
    open(INFILE, "<$seqFileName") or fatal("Could not open infile $seqFileName for reading.");
    $infh = *INFILE;
}

# open outfile or STDOUT
if(!defined($outFileName)) {
    $DEBUG && debug("Writing to STDOUT."); # default behaviour
} else {
    local *OUTFILE;
    open(OUTFILE, ">$outFileName") or fatal("Could not open $outFileName for writing.");
    $outfh = *OUTFILE;
}

# read in fasta sequences
$DEBUG && debug("Reading fasta sequences...");

my @seq;
my @seqName;
getFastaSeq($infh, \@seq, \@seqName);

# apply re to sequences, outputing results as we procede
print $outfh "Pattern $re matches\nMatch\tPos\tSequence name\n";

my $nr=0;
foreach $seqName (@seqName) {
    # match in forward seq
    while ($seq[$nr] =~ /$re/g) {
	# handle the returned matches
	$position = length($`) ; # length of prematch is position
	print $outfh "$&\t$position\t$seqName\n";
    }

    $nr++;
}

# end of main

sub debug {
    my $msg = shift;

    print STDERR "DEBUG: $msg\n";
}

sub fatal {
    my $msg = shift;

    print STDERR "FATAL: $msg\n";
    exit;
}

sub warning {
    my $msg = shift;

    print STDERR "WARNING: $msg\n";
}

sub usage {
    my $msg = shift;
    
    !defined($msg) && ($msg = "");
    $msg ne ""  && ($msg .= "\n");

    $msg .= "Usage: fastare [--infile <filename>][--outfile <filename>][--regexp <perl re>][--case-insensitive]\n";

    print STDERR $msg;
    exit;
}

sub getFastaSeq {
    my $fh = shift;
    
#  print "DEBUG: getFastaSeq tries to open $fastaFileName.\n";
#    open(FASTAFILE,"<$fastaFileName") || die "Sequence fasta input file $fastaFileName open failed.\n";
    
    my $seq = shift; # ref
    my $name = shift; # ref
    
    # First, get the sequences
    my $nFASTA=0;
    while(<$fh>) {
	chomp;
	if(/^\>/) {
	    # On fasta description line
	    $nFASTA++;
	    ($$name[$nFASTA-1])=/^\>(.+?)\s*$/; # skip trailing whitespaces
	    $$seq[$nFASTA-1]="";
	} else {
	    # Get all genomic sequence chars into that $seq string..
	    $$seq[$nFASTA-1].=$_;
	}
    }
    
    # done processing fasta sequence file - no close in case $fh == *STDIN
}
