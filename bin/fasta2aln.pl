#!/usr/bin/perl -w

use Bio::AlignIO;

my $inputfilename = $ARGV[0];
my $outputfilename = $ARGV[1];

my $in  = Bio::AlignIO->new(-file => $inputfilename ,
			 '-format' => 'fasta');
my $out = Bio::AlignIO->new(-file => ">$outputfilename" ,
			 '-format' => 'clustalw');

while ( my $aln = $in->next_aln() ) {
    $out->write_aln($aln);
}

