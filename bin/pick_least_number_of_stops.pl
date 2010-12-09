#!/usr/bin/perl -w

my $DEBUG = 1;

my $only_phrap_stopcount = $ARGV[0]; 
my $real_phrap_stopcount = $ARGV[1];
my $real_phrap_contignr_cutoff = $ARGV[2];

$DEBUG && print "running pick_least_number_of_stops with $only_phrap_stopcount $real_phrap_stopcount $real_phrap_contignr_cutoff.\n"; 

open ONLY_PHRAP, $only_phrap_stopcount;
open REAL_PHRAP, $real_phrap_stopcount;

open PHRAP_OUT,">$only_phrap_stopcount.better_junk_in_phrap";
open REAL_OUT,">$real_phrap_stopcount.better_junk_after_realign";

my %only_phrap;

while(my $r = <ONLY_PHRAP>) {

    chomp $r;
    $r =~ s/^\s+//;
    my @r = split(/\s+/,$r);

    my ($contig) = ($r[1]=~/([cC]ontig\d+)_/);
    my $stopcount = $r[0]; 
#    $DEBUG && print "only_p $contig, $stopcount\n";
    $only_phrap{$contig} = $stopcount;
}

my %real_phrap;

while(my $r = <REAL_PHRAP>) {

    chomp $r;
    $r =~ s/^\s+//;
    my @r = split(/\s+/,$r);

    my ($contig) = ($r[1]=~/([cC]ontig\d+)_/);
    my $stopcount = $r[0];
#    $DEBUG && print "real_p $contig, $stopcount\n";
    $real_phrap{$contig} = $stopcount;
}



my @phrap_only_contigs = keys %only_phrap;
my @real_phrap_contigs = keys %real_phrap;

my %uniq;
my $c; 
my @nonuniq = (@phrap_only_contigs,@real_phrap_contigs);
while ($c = shift @nonuniq) {
    $uniq{$c} = 1;
}
my @phrap_contigs = keys %uniq;

my %choose_contig_from;

while ($contig = shift @phrap_contigs) {    

    if(!defined($real_phrap{$contig})) {
	$real_phrap{$contig} = 1000000;
    }

    if(!defined($only_phrap{$contig})) {
	$only_phrap{$contig} = 1000000;
    }

    $DEBUG && print "DEBUG: $contig, ",$only_phrap{$contig}," ",$real_phrap{$contig},"\n";

    if($only_phrap{$contig} == 0) {       
	$choose_contig_from{$contig} = "phrap";	
    }elsif($only_phrap{$contig} < $real_phrap{$contig}) {
	$choose_contig_from{$contig} = "phrap";	
    } elsif ($only_phrap{$contig} > $real_phrap{$contig}) {
	$choose_contig_from{$contig} = "real";
    } elsif ($only_phrap{$contig} == $real_phrap{$contig}) {
	my ($contig_nr) = ($contig =~ /ontig(\d+)/);
	if($contig_nr <= $real_phrap_contignr_cutoff ) {
	    $choose_contig_from{$contig} = "phrap";
	} else {
	    $choose_contig_from{$contig} = "real";
	}
    } else {
	print "Someone has made a mess.\n";
	exit;
    }

    if($choose_contig_from{$contig} eq "phrap") {

	if($only_phrap{$contig} < 2) {
	    print PHRAP_OUT $contig."\n";
	} else {
	    $DEUG && print "Ignoring $contig. Best stopcount ",$only_phrap{$contig}," for only phrap.\n";
	}
    } else {
	if($real_phrap{$contig} < 2) {
	    print REAL_OUT $contig."\n";
	} else {
	    $DEUG && print "Ignoring $contig. Best stopcount ",$real_phrap{$contig}," for realigned contig.\n";
	}

    }
}
