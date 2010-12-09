#!/usr/bin/perl -w 

# simplistic, without any multifield feature..

my $field = 0;

while (my $arg = shift @ARGV) {
    if($arg =~ /-f/) {
	$field = shift @ARGV;
	if( $field < 1 ) {
	    print "Error. No field $field possible.";
	    exit;
	}
	$field--;
    }
}

while (<STDIN>) {
    my @row = split /\s+/;   
    if ($field < @row) {
	print $row[$field]."\n";
    }
 }
