#!/usr/bin/perl -w

my $DEBUG = 0;

$contignamefile = $ARGV[0];
$readlistfile = $ARGV[1];

open READLIST, "<$readlistfile";
open CONTIGNAMES, "<$contignamefile";

while (<CONTIGNAMES>) {
    chomp;
    push @contiglist, $_;
    $DEBUG && print "DEBUG: Got $_ for contiglist.\n";
}
close CONTIGNAMES;

my $in_printing_contig = 0;	
while($r = <READLIST>) {
   chomp $r;
    if( $r =~ /^CO (Contig\d+)/ ) {
	my $current_contig = $1;
	$DEBUG && print "DEBUG: found contig $current_contig\n";

	if ($in_printing_contig == 1) {
		$in_printing_contig = 0;
		# could remove from contiglist..
	}

        foreach my $contig ( @contiglist ) {
		if ( $current_contig eq $contig ) {

#			print $contig, "\n";
			open OUT, ">$contig.reads";
			$in_printing_contig = 1;
		}
	}

    } elsif ( $r =~ /^AF (\S+)/ ) {
	if ( $in_printing_contig ) {
		print OUT $1, "\n";
	}
    }

}

