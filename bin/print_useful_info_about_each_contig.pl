#!/usr/bin/perl -w

# first, among other things,
# count_fractions_in_ace_file.pl < se31reads.fasta.ace | sort -k2,2nr -k5,5nr
# then, run this

my $DEBUG = 0;
my $WARNING = 1;

my $DOMINANT_CUTOFF_DEFAULT = 4;

my $SE_UKS_COMPAT = 1;

my $dominant_cutoff = $DOMINANT_CUTOFF_DEFAULT;

my $dominance_exclusion_gbk_ids_file= "";
my $exclude = 0;
my $blasthitgff_file = "";
my $isolate_sort = 0;
my $dominants_only = 0;

my $typesort_file = "";
my $fractioncount_file = "";

my $co_dominants=0;
my $se_dominants=0;

# CLI

while (my $arg = shift @ARGV) {
    if($arg eq "--debug") {
	$DEBUG = 1;
    } elsif($arg eq "--typesort") {
	$typesort_file = shift @ARGV;
    } elsif($arg eq "--dominants_only") {
	$dominants_only = 1;
    } elsif($arg eq "--fractioncount") {
	$fractioncount_file = shift @ARGV;
    } elsif($arg eq "--dominant_cutoff") {
	$dominant_cutoff = shift @ARGV;
    } elsif($arg eq "--blasthitgff") { 
	$blasthitgff_file = shift @ARGV;
    } elsif($arg eq "--exclude_from_dominance") { # required for printing!
	$dominance_exclusion_gbk_ids_file = shift @ARGV;
#    } elsif($arg eq "--tree") {
#	$tree_file = shift @ARGV;
    } elsif ($arg eq "--no-warn") {
	$WARNING = 0;
    } else {
	print STDERR "WARNING: Unrecognised command line parameter $arg. Continuing anyway.\n";
    }
}

if( $dominance_exclusion_gbk_ids_file ne "" ) {
    if( $blasthitgff_file eq "") {
	$WARNING && print "No blasthits file given so exclusion of hit ids from dominance computation will not be made.\n";    
	$exclude = 0;
    } else {
	$exclude = 1;
	open(EXCLUDE, "<$dominance_exclusion_gbk_ids_file");
	while (<EXCLUDE>) {
	    chomp;
	    if( $_ ne "" ) {
		print STDERR "Exclude ",$_, "\n";
		push @exclude_list, $_;
	    }
	}
	close EXCLUDE;
    }
}

if($typesort_file eq "" || $fractioncount_file eq "") {
    usage();
}

open(TYPESORT, "<$typesort_file");
open(FRACTCOUNT, "<$fractioncount_file");
# open(TREE, "<$tree_file");
open(BLASTGFF, "<$blasthitgff_file");

my %typesort;

while(<TYPESORT>) {
    chomp;
    my (@row) = split /\s+/,$_;

    if( exists($typesort{$row[0]}) ) {
	$typesort{$row[0]} .= "_".$row[1]."_".$row[2];
    } else {
	$typesort{$row[0]} = $row[1]."_".$row[2];
    }
}

my %blasthit;
my %lowestexpect;

while( my $row = <BLASTGFF> ) {

    my @col = split /\t+/, $row; 

    my $contigname = $col[0];
    my $querystart = $col[3]; 
    my $queryend = $col[4];
    my $colelements = scalar(@col)-1;
    my $comment = join "\t", @col[8..$colelements];

    my @comment = split /;+/, $comment;
    
 #   print STDERR $comment[0]."\n";
    
    my $subjectname;
    if ($comment[0] =~ m/b\|([\w\d\.]+)\|/) {	
	$subjectname = $1;
    } else {
	($subjectname) = ($comment[0] =~ m/([\w\d\_\.]+)\|mRNA\|/);
    }

    my ($expect) = ($comment[1] =~ m/Expect = ([\deE\.\-]+)/);
    my ($querycoverage) = ($comment[1] =~ m/Querycoverage = ([\d\.]+)/);
    my ($idfreq) = ($comment[1] =~ m/Ids = \d+\/\d+ \(([\d\.]+)\)/);

    if(!defined($lowestexpect{$contigname}) || ($expect < $lowestexpect{$contigname})) { 
	$blasthit{$contigname} = "_".$subjectname."_$querystart+".$queryend."_$querycoverage"."_".$idfreq;
	$lowestexpect{$contigname} = $expect;
    }

#    $DEBUG && print STDERR $contigname." ".$blasthit{$contigname}." ".$lowestexpect{$contigname}."\n";    
}

my %exclude;

if($exclude == 1) {
    $DEBUG && print STDERR "checking excludes for ", scalar(keys %blasthit), " hit contigs.\n";
    
    # check exclusion id list fo to 
    foreach my $contig (keys %blasthit) { 

	$exclude{$contig} = 0;

	foreach $exclude_id (@exclude_list) {
	    if ($blasthit{$contig} =~ m/\_$exclude_id\_/) {
		$DEBUG && print STDERR "Excluding hit for $contig to excluded id $exclude_id.\n";
		$exclude{$contig} = 1;
	    }
	}
    }
}

my %isolate_fraction;
my %dominant;

while(<FRACTCOUNT>) {

    chomp;
    my @row = split /\s+/, $_;

    my $isolate = $row[2];
    my $contig = $row[0];
#    my $isolate_fraction = $row[5];
    my $isolate_fraction = $row[7];
    ${$isolate_fraction{$isolate}}->{$contig} = $isolate_fraction;

#    $DEBUG && print "DEBUG: $isolate $contig $isolate_fraction\n";
}

my %fractcount;

foreach my $isolate (keys %isolate_fraction) {

    my (@contigs) = keys %${$isolate_fraction{$isolate}};

    if ( $exclude == 1 ) {
	my @include_contigs;
	foreach my $contig (@contigs) {
	    if ( defined( $exclude{$contig} ) ) {
		if( $exclude{$contig} == 1) {
		    $DEBUG && print STDERR "Excluded contig $contig (hit ".$blasthit{$contig}.").\n";
		} elsif ($exclude{$contig} == 0) {
		    push @include_contigs, $contig;
		}
	    } else {
		# if contig didn't get a blast hit, mark it non-excluded =)
		$exclude{$contig} = 0; 
		# AND includet, you nipwit
		push @include_contigs, $contig;
	    }
	}
	@contigs = @include_contigs;
    }
    
    @contigs = sort { ${$isolate_fraction{$isolate}}->{$b} <=> ${$isolate_fraction{$isolate}}->{$a}} @contigs;

    my $maxcount = ( @contigs < $dominant_cutoff ) ? @contigs : $dominant_cutoff;
    
    $DEBUG && print "$maxcount dominants to process for $isolate (out of ",scalar(@contigs)," contigs).\n";

#    $DEBUG && print "DEBUG: $isolate maxcount $maxcount\n";
    
    for( my $i = 0; $i < $maxcount ; $i++ ) {
	$dominant{$contigs[$i]} = 1;
	if( exists($fractcount{$contigs[$i]}) ) {
	    # this is a (sub)dominant fraction
	    $fractcount{$contigs[$i]} .= "_".$isolate.".".($i+1)."_".${$isolate_fraction{$isolate}}->{$contigs[$i]};
#	    $DEBUG && print "DEBUG: ".$contigs[$i]." ".$fractcount{$contigs[$i]}."\n";
	} else {
	    $fractcount{$contigs[$i]} = $isolate.".".($i+1)."_".${$isolate_fraction{$isolate}}->{$contigs[$i]};
	    
#	    $DEBUG && print "DEBUG: ".$contigs[$i]." ".$fractcount{$contigs[$i]}."\n";
	}
    }
}

my ( @contigs ) = keys %exclude;
my %isolatetype_dominants;

foreach my $contig ( @contigs ) {

    if($dominants_only && !$dominant{$contig}) {
	next;
    }

    print "$contig\t";

    if(defined($blasthit{$contig})) {
	print "\t".$blasthit{$contig};
    } else {
	$DEBUG && print STDERR "No hits for $contig.\n";
    }

    if(exists($fractcount{$contig})) {
	print "\t",$fractcount{$contig};

       # ..and count this for the co/se random draw statistics..
	$fractcount = $fractcount{$contig};

	while ($fractcount =~ m/([a-zA-Z]+)(\d+)/g ) {
	    
	    my $isolatetype = $1;
#	    my $isolate = $1.$2;
	    

	    if ($SE_UKS_COMPAT) {

		# old, presumably redundant check for error in naming convention
		
		if( $isolatetype eq "cp" ) {
		    $isolatetype = "co";
		}

		# group uks & se for final tally?		
		if ( $isolatetype eq "uks" ) {
		    $isolatetype = "se";
		}
	    }
	    
	    $isolatetype_dominants{$isolatetype}++;
	}
	
	
# 	while ( $fractcount =~ m/co/g )  {
# 	    $co_dominants++;
# 	}

# 	while ( $fractcount =~ m/cp/g ) {
# 	    $co_dominants++;
# 	}
    
# 	while ( $fractcount =~ m/uks/g ) {
# 	    $se_dominants++;
# 	}
    
# 	while ( $fractcount =~ m/se/g ) {
# 	    $se_dominants++;
# 	}

    } 
    if(exists($typesort{$contig})) {
	print "\t".$typesort{$contig};
    }

    print "\n";
}

print STDERR "\nStats:";
foreach my $isolatetype ( keys %isolatetype_dominants ) {
    print STDERR " ",uc($isolatetype)," dominants: ", $isolatetype_dominants{$isolatetype};
}
print STDERR ".\n\n";

# CO dominants: $co_dominants SE dominants: $se_dominants.\n\n";
# nb: case may now change in the event of lc in.. setting to "uc" for old times sake..

sub usage {
    print STDERR "Usage: print_useful_info_about_each_contig.pl\n";
    print STDERR "\t<--typesort file>\n";
    print STDERR "\t<--fractioncount file>\n";
    
    print STDERR "\t[--dominant_cutoff integer]\t$DOMINANT_CUTOFF_DEFAULT\n";
    print STDERR "\t[--dominants_only]\t\tIgnore non-dominant contigs for final count.\n";

    print STDERR "\t[--blasthitgff file]\n";
    print STDERR "\t[--exclude_from_dominance file_with_gbkids]\n";

    exit;
}
