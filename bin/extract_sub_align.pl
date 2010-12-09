#!/usr/bin/perl -w

use Bio::AlignIO;

my @clargs = @ARGV;
my $DEBUG = 0;
my $WARNING = 1;

my $outfile = "";
my $infile = "";
my $sub_aln_start_col = -1;
my $sub_aln_end_col = -1;
my $rows_string = "";
my $inc_rows = 0;
my $del_rows = 0;

while (my $arg = shift @clargs) {

    if($arg eq "--out") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no filename given to --out flag.\n";
	} else {
	    $outfile = $arg;
	}
    } elsif ($arg eq "--in") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no filename given to --in flag.\n";
	} else {
	    $infile = $arg;
	}
    } elsif($arg eq "--start") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no argument given to --start flag.\n";
	} else {
	    $sub_aln_start_col = $arg;
	}
    } elsif($arg eq "--end") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no argument given to --end flag.\n";
	} else {
	    $sub_aln_end_col = $arg;
	}
    } elsif($arg eq "--delrows") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no rows to delete given to --delrows flag.\n";
	} else {
	    $rows_string = $arg;
	    $del_rows = 1;
	}
    } elsif($arg eq "--rows") {
	$arg = shift @clargs;
	
	if($arg eq "") {
	    $WARNING && print STDERR "WARNING: no rows to include given to --rows flag.\n";
	} else {
	    $rows_string = $arg;
	    $inc_rows = 1;
	}
    } elsif($arg eq "--debug") {
	$DEBUG = 1;
	$WARNING && print STDERR "DEBUG mode enabled.\n";
    } elsif($arg eq "--silent") {
	$WARNING = 0;
    } else {
	$WARNING && print STDERR "WARNING: Unrecognised command line parameter $arg.\n";
	usage();
    }

}

if($del_rows && $inc_rows) {
    $WARNING && print STDERR "FATAL: cannot both delete and include rows. Please choose one mode.\n";
}

if( $outfile eq "" or $infile eq "" or $sub_aln_end_col == -1 or $sub_aln_start_col == -1 ) {
    $WARNING && print STDERR "FATAL: Missing parameter.\n";
    usage();
}

my @rows = ();
if( $rows_string ne "" ) {
    my @rows_elem = split(/,/, $rows_string);
    while (my $element = shift @rows_elem) {
	if ( ($e_start, $e_end) = ( $element =~ m/(\d+)-(\d+)/ )) {
	    push @rows, $e_start..$e_end;
	} else {
	    push @rows, $element;
	}
	@rows = sort { $a <=> $b } @rows;
    }

    $DEBUG && print STDERR "rows: ", join(" ", @rows),"\n";
}

$ain = Bio::AlignIO->new('-file'=>$infile, '-format'=>'clustalw');
$aout = Bio::AlignIO->new('-file'=>">$outfile", '-format'=>'clustalw');

my $aln = $ain->next_aln();

if ($del_rows) {
    while ( my $r = shift @rows ) {
	my $tmp_seq = get_seq_by_pos($r);
	$aln->remove_seq($tmp_seq);
    }
} elsif($inc_rows) {   
#    my $tmp_seq = $aln->get_seq_by_pos($rows[0]);

#    my $start_loc = $tmp_seq->location_from_column( $sub_aln_start_col );
#    my $end_loc = $tmp_seq->location_from_column( $sub_aln_end_col );
    $aln = $aln->select_noncont( @rows );
#    $sub_aln_start_col = $aln->column_from_residue_number($tmp_seq->id, $start_loc->start );
#    $sub_aln_end_col = $aln->column_from_residue_number($tmp_seq->id, $end_loc->start );
}

# to get the slice pos right, save col pos coord in one of the included, and recapture column..
my $slice = $aln->slice($sub_aln_start_col, $sub_aln_end_col);

$aout->write_aln($slice);

sub usage {
    print "USAGE:\textract_sub_align <--in clustalw.aln>";
    print "\t\t<--out clustalw.aln>\n";
    print "\t\t<--start start.col>\n";
    print "\t\t<--end end.col>\n";
    print "\t\t[--delrows rows,to-remove,from-out,aln | --rows rows,to-include,from-out,aln]\n";

    exit(1);
}
