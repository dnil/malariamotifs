#!/usr/bin/perl -w

use strict;

my $DEBUG = 0;

my $acefilename = $ARGV[0];
my $cid = $ARGV[1];

if (@ARGV < 2 ){
    print STDERR "USAGE: draw_contig.pl acefile contigid > page.html\n";
    exit;
}

# globals
my %readseq;
my %readstart;
my %readlen;
my %matename;
my $contig_seq = "";

my $mostnegative = 1;

sub read_ace;

read_ace($acefilename, $cid);

## Draw! ##
my $initspaces = 0;
# end of interesting contig. print and exit.

my $contig_length = length($contig_seq);

if ($mostnegative < 1) {
    $initspaces= -$mostnegative +1; # 0 possibe? if so, add another space..
#    if ($mostnegative == 0) {
#	print "OOPS: mostnegative is 0. fix!\n";
#    }
}

my $view_length = $contig_length + $initspaces ;

print "<html><head><title>$acefilename - $cid</title>";
print "<STYLE type=\"text/css\"><!--\n";
print ".basicread { color: black; }\n";
print ".highlightread { color: red; background-color: maroon}\n";
print ".mateless {color: purple; }\n";
print ".consensus {background-color: silver; }\n";
print ".ruler {color: olive;}\n";
print ".tooltip {visibility: hidden; position: absolute; border-color: black; border-width: 1; background-color: yellow; color: black}";
print ".spacer {}\n";
print "-->\n</STYLE>";

print "<meta http-equiv=\"Content-Script-Type\" content=\"text/javascript\">\n";
#print "<script>";
#print "";
#print "</script>";

print "</head><body>\n";

print "<div class=\"tooltip\" id=\"tooltip\"></div>"; # a generic tooltip, to move and make visible as needed.

#make ruler line 
my $rulerline = "<tt class=\"ruler\">";

$rulerline .= "&nbsp;"x$initspaces;
for (my $i=100; $i<$contig_length; $i += 100) {
    my @ruler_val = split(/ */, "$i");

    $rulerline .= "."x(100 -1 - @ruler_val );
    $rulerline .= "$i|";
}
$rulerline .= "</tt>\n";

# and print once
print $rulerline;

print "<tt id=\"contig\" class=\"consensus\">";
print "&nbsp;"x$initspaces;
print $contig_seq;
print "</tt>\n";

my $count = 0;

# sort reads according to readstart

my @reads = keys %readstart; 

my @sorted_reads = sort { $readstart{$a} <=> $readstart{$b} } @reads;

foreach my $read (@sorted_reads) {
    
    #give a ruler every N reads
    $count++;
    if ($count % 50 == 0) {
	print $rulerline;
    }

    my $read_initspaces = $initspaces + $readstart{$read}-1;
#    if ($readstart{$read}<1) {
#	$read_initspaces -= 1;
#    }
    print "<tt class=\"spacer\">";
    print "&nbsp;"x$read_initspaces;
    print "</tt>";

    my $html_id = $read;
    $html_id =~ s/\./_/g;

    if(defined($matename{$read})) {
	my $mate_id = $matename{$read};
	$mate_id =~ s/\./_/g;

	my $mouseover = "this.className=&quot;highlightread&quot;; document.getElementById('$mate_id').className=&quot;highlightread&quot;; document.getElementByID('tooltip').innerHTML='<p>$html_id ($mate_id)<\/p>'; document.getElementByID('tooltip').left=event.pageX+&quot;px&quot;; document.getElementByID('tooltip').top=event.pageY+&quot;px&quot;; document.getElementByID('tooltip').style.visibility=&quot;visible&quot;;";

	my $onclick = "document.getElementById('$mate_id').scrollIntoView(true); document.getElementByID('tooltip').style.visibility='hidden';";
	my $mouseout = "this.className=&quot;basicread&quot;; document.getElementById('$mate_id').className=&quot;basicread&quot;; document.getElementByID('tooltip').style.visibility=&quot;hidden&quot;;";
	print "<tt id=\"$html_id\" onmouseover=\"$mouseover\" onmouseout=\"$mouseout\" onclick=\"$onclick\">"; # class=\"basicread}\"
    } else {
	print "<tt id=\"$html_id\" class=\"mateless\">";
    }       
    print $readseq{$read},"[$read]";
    print "</tt>\n";

}

print "</body></html>\n";

sub read_ace {
    my $acefilename = shift;
    my $cid = shift;

    my $correct_contig=0;
    my $in_contig_seq;

#    my $current_contig;
    my $current_contig_len;
    my $current_contig_reads;
        
    open ACE, "<$acefilename";
        
    my $in_read = 0;
    my $current_read = "";

    while(my $acer= <ACE>) {
	chomp $acer;
	
	if($correct_contig ==1) {		
	    if($in_contig_seq == 1) {
		if($acer eq "") {
		$in_contig_seq = 0;
		$DEBUG && print "done parsing seq\n";
	    } else {
		$DEBUG && print "Adding $acer to seq..\n";
		$contig_seq .= $acer;
	    }
	}
	
	if($in_read == 1) {
	    if( $acer eq "" ) {
		# done
		$in_read = 0;
		$readlen{$current_read}=length( $readseq{$current_read} );		
#		if($readlen{$current_read}+$readstart{current_read} > $current_contig_len) {
#		    $longest_read_overflow = $readlen{$current_read}+$readstart{current_read} - $current_contig_len; # assuming no individual read will be longer than the contig itself
#		}

		$DEBUG && print "done parsing $current_read\n";
	    } else {
		$readseq{$current_read} .= $acer;
		$DEBUG && print "Adding $acer to read seq for $current_read.\n";
	    }
	}
	
	if($acer =~ m/^AF (\S+) [UC] ([-\d]+)/) {
	    $current_read = $1;
	    
	    $readstart{$current_read} = $2;
	    $DEBUG && print STDERR "Found read $current_read at pos $2..\n";
	    if ( $readstart{$current_read} < $mostnegative) {
		$mostnegative = $readstart{$current_read};
	    }
	    
	    # does it have any mate yet?
	    # putative mate name

	    $current_read =~ m/^([frFR])([^.]+)\.?(\d+)\.(.+)/;

	    my $orientation = $1;
	    my $read_base = $2;	   
	    my $clone_nr = $3;
	    my $remainder = $4;

	    my $rev_orientation ="";
	    if ($orientation eq "f") {
		$rev_orientation = "r";
	    } elsif ($orientation eq "F") {
		$rev_orientation = "R";
	    }elsif ($orientation eq "r") {
		$rev_orientation = "f";
	    }elsif ($orientation eq "R") {
		$rev_orientation = "F";
	    }
	    my $putative_matename = $rev_orientation.$read_base.".".$clone_nr.".".$remainder;
	    $DEBUG && print STDERR "putative mate for $current_read is $putative_matename\n";    
	    if( defined($readstart{$putative_matename}) ) {
		$DEBUG && print STDERR "found mate!\n";
		$matename{$current_read} = $putative_matename;
		$matename{$putative_matename} = $current_read;
	    }
	}

	if ($acer =~ m/^RD (\S+) (\d+) \d+ \d+/) {
	    $current_read = $1;
#	my $current_read_len = $2;
	    
	    $readseq{$current_read} = "";
	    $in_read = 1;
	}

	if ($acer =~ m/^CO (Contig\d+)/ or $acer =~ /^WA/) {
	    $DEBUG && print STDERR "done parsing contig.\n";
	    
	    return;

	}
    }
    
    if($acer =~ m/^CO (Contig\d+) (\d+) (\d+)/) {
	my $current_contig = $1;
	$current_contig_len = $2;
	$current_contig_reads = $3;

	
	if($current_contig eq $cid) {
	    $DEBUG && print STDERR "Found $current_contig...\n";
	    $in_contig_seq = 1;
	    $correct_contig = 1;
	} else {
	    $correct_contig = 0;
	}
    }
    }

}

