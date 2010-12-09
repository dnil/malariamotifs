#!/usr/bin/perl -w

my $optima = -1;
my $scoretype = "S";

while (my $arg = shift @ARGV) {
    if ($arg eq "--min") {
	$optima = -1;
    } elsif ($arg eq "--max") {
	$optima = 1;
    } elsif ($arg eq "--field" ) {
	$arg = shift @ARGV;
	$score_type = $arg;
    } else {
	print STDERR "Unknown option \"$arg\".\n";
	usage();
	exit;
    }
}

my %domn;
my %ids; 
while ($row = <STDIN>) {
    chomp $row; 
    ($score) = ($row =~ /$scoretype = ([-+\d\.]+)/);
    (@isolates) = ($row =~ /_(\w{2,3}\d{1,2}\.\d)_/g);
    while ($i = shift @isolates) {
	my ($iso, $dom) = ($i =~/(\w{2,3}\d{1,2})\.(\d)/);

	if( !defined($ids{$iso}) || ($optima < 0 && $score < $ids{$iso}) || ($optima > 0 && $score > $ids{$iso}) ) {
	    $ids{$iso} = $score;
	    $domn{$iso} = $dom;
	}
    }
} 

@keys = sort keys %ids; 
while ($k = shift @keys) {
    print $k, "\t", $ids{$k},"\t", $domn{$k},"\n";
}


sub usage {
    print STDERR "Usage: query.pl < <STDIN> > <STDOUT> [<--min>|--max] [--field <S>|M|SE]\n";
}
