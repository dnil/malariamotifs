#!/bin/bash
#
# process_reads.sh read_file_name task [matrix]
#
# task, the second parameter, is one of 
#
#      all        -  all of the below
#      phred      -  basecall reads (requires phred, www.phrap.org)
#      trim       -  quality trim (requires EMBOSS, www.emboss.org)
#      phrap      -  assembly/clustering (requires phrap, www.phrap.org)
#      translate  -  conceptual translation (requires EMBOSS)
#      realign    -  realigner (requires realigner, binary shipped else enquire from bjorn.andersson@ki.se)
#                    (also translates realigned reads, req. EMBOSS)
#      flip       -  reorient nucleotide sequences based on best peptide trans using known motifs and top codon count (req EMBOSS)
#      blast      -  attempt to annotate and recognise known sequences (req NCBI blast and a formatted database, preferenially from the gene of interest)    
#                    also instrumental in filtering out possible contaminants and spearately treating e.g. possible known untranslated/unexported genes.
#      tree       -  tree construction (req clustal, www.clustal.org)
#      motifs     -  (req R, www.r-project.org and NCBI blast, EMBOSS for PSSM-export and finding those motifs in other sequences) 
#                    Also run section_alignment w/o options for usage. See also section alignment v2 and rescore using R for more details.
#      cleanup    -  attempts to clean up temporary files etc, but as the concept of temporary is somewhat depending on the task at hand, this is rather permissive.
#      
# The aa-similarity matrix, optional third argument, can be either blosum or gonnet.
#
# There are also a couple of other tasks, lmo and stats, that are helpful for development etc. Can also be instructive for how to classify the found motifs in new sequences etc. 
#
# The required binaries are assumed to be on the PATH, so please set it accordingly.
# 
# I don't like to apologize about code, but this was an unusually evolutionary project and written way back when. I guess it's a sign of progress that I 
# at least now know many, many ways to make this better now. Actually probably already could think of many back then, but time in PhD sideprojects
# is a complicated resource to manage. If someone could get me a couple of weeks I could get this up to a distributable shape.. 
#
# Daniel Nilsson, 2004-2005. 
#

READFILE=$1;

READFILE_BASENAME=`basename $READFILE .fasta`

BASE_DIR="/home/daniel/pf"

export PATH=$PATH":/home/daniel/src/phrap:/home/daniel/src/blast-2.2.24/bin:/home/daniel/src/phred.071220b"
export PHRED_PARAMETER_FILE=/home/daniel/src/phred.071220b/phredpar.dat

BIN_DIR=$BASE_DIR"/malariamotifs/bin"
DATA_DIR=$BASE_DIR"/malariamotifs/data"
PRIMER_PATTERN_FILE=$DATA_DIR"/primer_patterns.txt"
DOMINANCE_EXCLUDE=$DATA_DIR"/var_common_gbk_ids"

BLASTDB="pf_gbk_plasmodb_040512_var_gene_sequences.fasta"
DOMINANT_CUTOFF=3;

# temporarily enabled.. running "most" only..
#READFILE_BASENAME=malaria.040526.screen.qp.seq.namefix.loq

#ACEFILE=malaria.040526.screen.qp.seq.namefix.loq.ace
ACEFILE=$READFILE_BASENAME".screen.qp.plq.seq.ace"

REAL_CONS_BASENAME=$READFILE_BASENAME".real.cons"
ACEBASE=`basename $ACEFILE .ace`

# big! 
#DB_NAME=$READFILE_BASENAME.par.clustalnames.nonsingletons.pep.aln.db.fasta
# "Small-"!
#DB_NAME=$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.aln.db.fasta
# definitions exist in the motif-section already..

# process_reads.sh read_file_basename phred
# process_reads.sh read_file_basename trim
# process_reads.sh read_file_basename phrap
# process_reads.sh read_file_basename translate
# process_reads.sh read_file_basename realign
# process_reads.sh read_file_basename flip
# process_reads.sh read_file_basename blast
# process_reads.sh read_file_basename tree [<blosum>|gonnet]
# process_reads.sh read_file_basename motifs [<blosum>|gonnet]
# process_reads.sh read_file_basename cleanup 

ARSFADIG_MOTIF=RSFADIGDI
#LA are note always detected for the AF-primers..
STATS=$READFILE_BASENAME".stats.log"

ENTRYPOINT=$2
# goto $ENTRYPOINT..
echo $ENTRYPOINT" for "$READFILE

case "$ENTRYPOINT" in ( phred )
# basecall w/ phred
phred -id ../chromat_dir -pd ../phd_dir -sa $READFILE_BASENAME.seq -qa $READFILE_BASENAME.seq.qual

;;
( "trim" )
# cross_match - vector mask
cross_match -minmatch 10 -minscore 20 -screen $READFILE_BASENAME.seq $DATA_DIR/pcriitopo_seq.txt
#cross_match -minmatch 10 -minscore 20 -screen $READFILE_BASENAME.seq $DATA_DIR/pcr2_1topo_seq.txt
#cross_match -minmatch 10 -minscore 20 -screen $READFILE_BASENAME.seq $DATA_DIR/pcr4topo_seq.txt

# quality trim and trim the masked vector

$BIN_DIR/pruneESTsForClustering_x_mod.pl $READFILE_BASENAME".seq.screen" $READFILE_BASENAME".seq.qual" $READFILE_BASENAME".screen.qp.seq" --prune-x --write-qual $READFILE_BASENAME".screen.qp.seq.qual" > $READFILE_BASENAME".screen.qp.log"

# loqual primerresidues

#$BIN_DIR/loqual_ambigous_primer_bases.pl $PRIMER_PATTERN_FILE $READFILE_BASENAME.screen.qp.seq $READFILE_BASENAME.screen.qp.seq.qual $READFILE_BASENAME.screen.qp.plq.seq.qual > $READFILE_BASENAME.plq.log
$BIN_DIR/loqual_entire_primer_region.pl $PRIMER_PATTERN_FILE $READFILE_BASENAME.screen.qp.seq $READFILE_BASENAME.screen.qp.seq.qual $READFILE_BASENAME.screen.qp.plq.seq.qual > $READFILE_BASENAME.plq.log

# loqual gives new qual file so we need a new name for the corresponding sequence file
ln -s $READFILE_BASENAME".screen.qp.seq" $READFILE_BASENAME".screen.qp.plq.seq"

;;
( "phrap" )
# run phrap 

phrap -new_ace -minmatch 20 -retain_duplicates -repeat_stringency 0.9 $READFILE_BASENAME".screen.qp.plq.seq"  > $READFILE_BASENAME".phrap_log"

#ACEFILE=$READFILE_BASENAME".screen.qp.plq.seq.ace"
#ACEBASE=$READFILE_BASENAME".screen.qp.plq.seq"
#grep -A3 WA $ACEFILE |tail +2 >> $STATS

# get singletons

;;
( "translate" )

grep Contig $ACEFILE | cut -f2,4 -d' ' | sort -k2n | perl -ne 'if(m/\ 1$/) { } else {print;}' |cut -f1 -d' ' > $READFILE_BASENAME.nonsingletons.list

#$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.nonsingletons.list $READFILE_BASENAME.contigs > $READFILE_BASENAME.nonsingletons.fasta
$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.nonsingletons.list $ACEBASE.contigs > $READFILE_BASENAME.nonsingletons.fasta
transeq -frame 6 -sequence $READFILE_BASENAME.nonsingletons.fasta -outseq $READFILE_BASENAME.nonsingletons.pep

perl -ne 'chomp; if (m/^>.+(Contig\d+_\d+)/) { print ">",$1,"\n"; } else { print $_,"\n"; }' < $READFILE_BASENAME.nonsingletons.pep > $READFILE_BASENAME.nonsingletons.clustalnames.pep
sed -e 's/\*//g;' < $READFILE_BASENAME.nonsingletons.clustalnames.pep > $READFILE_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed

fuzzpro -pattern $ARSFADIG_MOTIF -pmismatch 3 -sequence $READFILE_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed -outfile $READFILE_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.fuzzpro

grep Sequence: $READFILE_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.fuzzpro | cut -f 3 -d' ' > $READFILE_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.list.unedited

cut -f1 -d'_' $READFILE_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.list.unedited |uniq -c |sort |grep \ 2\  > $READFILE_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes.1
$BIN_DIR/wscut.pl -f 3 < $READFILE_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes.1 > $READFILE_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes
grep -A9 -f $READFILE_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes $READFILE_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.fuzzpro |perl -e '$best = 4; $pair = 0; $last_match = ""; $current_match = ""; while (<STDIN>) { if (m/Sequence:\s+(Contig\d+_\d+)\s+/) { $last_match = $current_match; $current_match= $1; $pair = (($pair==0) ? 1 : 0); $DEBUG && print "DEBUG $pair $current_match\n"} elsif (m/End/) { $next_row = <STDIN>; ($mismatches_this_row) = ($next_row =~ m/^\s+\d+\s+\d+\s+([\.0-9]+)/); $mismatches_this_row == "." && ($mismatches_this_row = 0); if($pair == 0) { if($mismatches_this_row < $best) { print $last_match,"\n"; } else { print $current_match,"\n"; } } else { $best = $mismatches_this_row; } } }' >  $READFILE_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes.list

grep -v -f $READFILE_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes.list $READFILE_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.list.unedited > $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list

# next is a slow grep step; get all translated ARSFADIG-positive frames into one file..
$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list $READFILE_BASENAME.nonsingletons.clustalnames.pep > $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep

$BIN_DIR/fasta_pep_re.pl --regexp \\* --infile $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep |cut -f3 |grep Contig |sort -k1.7n |uniq -c |sort -n > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort

# pick NON-stop codon containing then
$BIN_DIR/wscut.pl -f 3 < $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort |sort > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list
grep -v -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list
$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep > $READFILE_BASENAME.clustalnames.nonsingletons.nostop.pep 

#echo -ne "\nNr of reads in stop codon containing ARSFADIG positive reading frames: "  >> $STATS
#cut -d_ -f1 $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.frameless
#grep -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.frameless $ACEFILE|cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

#echo -ne "\nNumber of reads in stop codon free ARSFADIG positive reading frames: " >> $STATS
#cut -d_ -f1 $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list.frameless
#grep -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list.frameless $ACEFILE|cut -f2,4 -d' '|awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

# realign, and keep both none and one stop frames w/ ARSFADIG, not precviously on the reoll

;;
( "realign" )

mkdir $READFILE_BASENAME.raws
mkdir $READFILE_BASENAME.real
mkdir $READFILE_BASENAME.real.fasta
mkdir $READFILE_BASENAME.aln.real
mkdir $READFILE_BASENAME.cons

cd $READFILE_BASENAME.raws

$BIN_DIR/ace2raw.pl --infile ../$ACEFILE
for file in *raw; 
  do 
    $BIN_DIR/Realigner -b 20 < $file > ../$READFILE_BASENAME.real/`basename $file .raw`.real
    $BIN_DIR/raw2multifasta.pl `basename $file .raw`.ids ../$READFILE_BASENAME.real/`basename $file .raw`.real > ../$READFILE_BASENAME.real.fasta/`basename $file .raw`.real.fasta
    $BIN_DIR/fasta2aln.pl ../$READFILE_BASENAME.real.fasta/`basename $file .raw`.real.fasta ../$READFILE_BASENAME.aln.real/`basename $file .raw`.real.aln
    $BIN_DIR/multifasta2cons.pl ../$READFILE_BASENAME.real.fasta/`basename $file .raw`.real.fasta ../$READFILE_BASENAME.cons/`basename $file .real.fasta`.real.cons.fasta
done

cd ..
cd $READFILE_BASENAME.cons
cat *real.cons.fasta > ../$READFILE_BASENAME.real.cons.fasta
cd ..

# translate and check for stop codons in the reals

perl -ne 'chomp; if (m/^>(\d+)/) { print ">Contig".$1,"\n"; } else { print $_,"\n"; }' < $REAL_CONS_BASENAME.fasta > $REAL_CONS_BASENAME.clustalnames.fasta
$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.nonsingletons.list $REAL_CONS_BASENAME.clustalnames.fasta > $REAL_CONS_BASENAME.nonsingletons.clustalnames.fasta
transeq -frame 6 -sequence $REAL_CONS_BASENAME.nonsingletons.clustalnames.fasta -outseq $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep

sed -e 's/\*//g;' < $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep > $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed

fuzzpro -pattern $ARSFADIG_MOTIF -pmismatch 3 -sequence $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed -outfile $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.fuzzpro

grep Sequence: $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.fuzzpro | cut -f 3 -d' ' > $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.list.unedited

cut -f1 -d'_' $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.list.unedited |uniq -c |sort |grep \ 2\  > $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes.1
$BIN_DIR/wscut.pl -f 3 < $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes.1 > $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes
grep -A9 -f $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.fuzzpro |perl -e '$best = 4; $pair = 0; $last_match = ""; $current_match = ""; while (<STDIN>) { if (m/Sequence:\s+(Contig\d+_\d+)\s+/) { $last_match = $current_match; $current_match= $1; $pair = (($pair==0) ? 1 : 0); $DEBUG && print "DEBUG $pair $current_match\n"} elsif (m/End/) { $next_row = <STDIN>; ($mismatches_this_row) = ($next_row =~ m/^\s+\d+\s+\d+\s+([\.0-9]+)/); $mismatches_this_row == "." && ($mismatches_this_row = 0); if($pair == 0) { if($mismatches_this_row < $best) { print $last_match,"\n"; } else { print $current_match,"\n"; } } else { $best = $mismatches_this_row; } } }' >  $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes.list

grep -v -f $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes.list $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.list.unedited > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list

# next is a slow grep step; get all translated ARSFADIG-positive frames into one file..
$BIN_DIR/fasta_header_grep.pl -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep >  $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep

$BIN_DIR/fasta_pep_re.pl --regexp \\* --infile $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep |cut -f3 |grep Contig |sort -k1.7n |uniq -c |sort -n > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort

$BIN_DIR/wscut.pl -f 3 < $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort |sort > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list

#echo -ne "\nNumber of reads in RE-AL stop codon containing ARSFADIG positive reading frames: " >> $STATS
#cut -d_ -f1 $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.frameless
#grep -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.frameless $ACEFILE |cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

# pick NON-stop codon containing then
grep -v -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nostop.list
$BIN_DIR/fasta_header_grep.pl -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nostop.list $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep > $REAL_CONS_BASENAME.clustalnames.nonsingletons.nostop.pep 

#echo -ne "\nNumber of reads in RE-AL stop codon free ARSFADIG positive contigs: " >> $STATS
cut -d_ -f1 $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nostop.list > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nostop.list.frameless
#grep -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nostop.list.frameless $ACEFILE |cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

#
# merge realigned set with original phrap set
#

# add 0s for the nostop realigned contigs to get a complete stopcount
cat $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.stopcount
perl -ne 'print "\t0 $_";' < $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nostop.list >> $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.stopcount

cat $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort > $READFILE_BASENAME.clustalnames.nonsingletons.pep.stopcount
perl -ne 'print "\t0 $_";' < $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list >> $READFILE_BASENAME.clustalnames.nonsingletons.pep.stopcount

# compare stop containing contigs to any realignment; if realigned has less stop codons, keep it
# perhaps this could be done once before ARSFADIG-matching as well, in order to save a few more reads?
# first 10read contig is a rather arbitrary cutoff for where the realigns with equal number of stops to the phrap contigs are considered better 
FIRST_10_READ_CONTIG=`grep ^CO $ACEFILE |cut -f2,4 -d' '|sed -e's/Contig//' |perl -ne 'chomp; @r = split(/\s+/); if ($r[1] == 10) { print $r[0]-1,"\n"; exit }'`
$BIN_DIR/pick_least_number_of_stops.pl $READFILE_BASENAME.clustalnames.nonsingletons.pep.stopcount $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.stopcount $FIRST_10_READ_CONTIG

# add an underscore to id without frame designation to desambiguate 1yx from 1yxz
perl -ne 'chomp; print $_."_\n";' < $READFILE_BASENAME.clustalnames.nonsingletons.pep.stopcount.better_junk_in_phrap > $READFILE_BASENAME.clustalnames.nonsingletons.pep.better_junk_in_phrap_
perl -ne 'chomp; print $_."_\n";' < $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.stopcount.better_junk_after_realign > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.better_junk_after_realign_

grep -i -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.better_junk_in_phrap_ $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.better_junk_in_phrap.list
grep -i -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.better_junk_after_realign_ $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.better_junk_after_realign.list

$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.better_junk_in_phrap.list $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep > $READFILE_BASENAME.clustalnames.nonsingletons.pep.better_junk_in_phrap.pep

$BIN_DIR/fasta_header_grep.pl -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.better_junk_after_realign.list $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.better_junk_after_realign.pep

# STATS preparatory.. probably not needed any longer?
#cut -d_ -f1 $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.frameless

#cut -f1 -d'_' $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.1
#$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.1 $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep > $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.with_stops_in_phrap_cons  

# cut -f1 -d'_' $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list.f1
# cut -f1 -d'_' $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.f1
# grep -v -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list.f1 $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.f1 > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.not_in_phrap_ok.list
# $BIN_DIR/fasta_header_grep.pl -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.not_in_phrap_ok.list $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.not_in_phrap_ok.pep

#cat $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.not_in_phrap_ok.pep $READFILE_BASENAME.clustalnames.nonsingletons.nostop.pep > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep
cat $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.better_junk_after_realign.pep $READFILE_BASENAME.clustalnames.nonsingletons.pep.better_junk_in_phrap.pep > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep

sed -e 's/\*//g;' < $READFILE_BASENAME.par.clustalnames.nonsingletons.pep > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed

fuzzpro -pattern KYSFADI -pmismatch 0 -sequence $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed -outfile  $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.KYSFADI.fuzzpro
grep Sequence: $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.KYSFADI.fuzzpro | cut -f 3 -d' ' >  $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.KYSFADI.list

fuzzpro -pattern RYSFADI -pmismatch 0 -sequence $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed -outfile $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.RYSFADI.fuzzpro
grep Sequence: $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.RYSFADI.fuzzpro | cut -f 3 -d' ' > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.RYSFADI.list

fuzzpro -pattern RYSFADL -pmismatch 0 -sequence $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed -outfile  $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.RYSFADL.fuzzpro
grep Sequence: $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.RYSFADL.fuzzpro | cut -f 3 -d' ' > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.RYSFADL.list

cat $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.RYSFADL.list $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.RYSFADI.list $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.KYSFADI.list > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.dbl2b.list

grep \> $READFILE_BASENAME.par.clustalnames.nonsingletons.pep |sed -e 's/>//;' > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.list

grep -v -f $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.dbl2b.list $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.list > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.list.without_dbl2b

#  $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.dbl2b.list > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.dbl2b.list

cp $READFILE_BASENAME.par.clustalnames.nonsingletons.pep $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.with_dbl2b

$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.list.without_dbl2b $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.with_dbl2b > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep

perl -ne 'chomp; if (m/^>.*(Contig\d+)/) { print ">",$1,"\n"; } else { print $_,"\n"; }' < $READFILE_BASENAME.par.clustalnames.nonsingletons.pep > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.clustalnames

# if more than one stop in best version, discard contig -- now done in pick_least_number_of_stops.pl
$BIN_DIR/fasta_pep_re.pl --regexp \\* --infile $READFILE_BASENAME.par.clustalnames.nonsingletons.pep |cut -f3 |grep Contig |sort -k1.7n |uniq -c |sort -n > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort

#echo -ne "\nNr of reads in stop codon containing ARSFADIG positive reading frames: "  >> $STATS
#cut -d_ -f1 $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.frameless

$BIN_DIR/wscut.pl -f 3 < $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort |sort > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort.list
# cut -d_ -f1 $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort.list > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort.list.frameless
#grep -f $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort.list.frameless $ACEFILE|cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

#
# flip the sequences with ARSFADIG on "-" strand (preparation for clustalw on the nt level)
#

# REINSERT unlucky contigs (like broken ARSFADIG but otw good, manually fixed fs, ...)

;;
( "flip" )

perl -ne 'chomp; if (m/^>.+(Contig\d+)/) { print ">",$1,"\n"; } else { print $_,"\n"; }' < $ACEBASE.contigs > $READFILE_BASENAME.clustalnames.contigs

grep \> $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep | sed -e 's/>//; s/_/ /;' > $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.frames

perl -ne 'chomp; @r = split /\s+/; if($r[1] >3) { print $r[0]."\n"; }' < $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.frames > $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.flipus

rm -f $READFILE_BASENAME.contigs.ARSFADIG.3mism.clustalnames.flipus.seq
for entry in `cat $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.flipus`; do $BIN_DIR/fasta_header_grep.pl $entry$ $READFILE_BASENAME.clustalnames.contigs >> $READFILE_BASENAME.contigs.ARSFADIG.3mism.clustalnames.flipus.seq ; done

revseq $READFILE_BASENAME.contigs.ARSFADIG.3mism.clustalnames.flipus.seq $READFILE_BASENAME.contigs.ARSFADIG.3mism.flipped

perl -ne 'chomp; @r = split /\s+/; if($r[1] <=3) { print $r[0]."\n"; }' <  $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.frames > $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.dont_flipus

rm -f $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.dont_flipus.seq
for entry in `cat $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.dont_flipus`; do $BIN_DIR/fasta_header_grep.pl $entry$ $READFILE_BASENAME.clustalnames.contigs >> $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.dont_flipus.seq ; done

cat $READFILE_BASENAME.contigs.ARSFADIG.3mism.flipped $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.dont_flipus.seq > $READFILE_BASENAME.contigs.ARSFADIG.3mism.flipped.seq

# flip the realingned sequences
grep \> $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep | sed -e 's/>//; s/_/ /;' > $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.frames

perl -ne 'chomp; @r = split /\s+/; if($r[1] >3) { print $r[0]."\n"; }' < $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.frames > $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.flipus

rm -f $READFILE_BASENAME.real.cons.fasta.flipus.seq
for entry in `cat $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.flipus`; do $BIN_DIR/fasta_header_grep.pl $entry$ $REAL_CONS_BASENAME.clustalnames.fasta >> $READFILE_BASENAME.real.cons.fasta.flipus.seq ; done

revseq $READFILE_BASENAME.real.cons.fasta.flipus.seq $READFILE_BASENAME.real.cons.fasta.flipped

perl -ne 'chomp; @r = split /\s+/; if($r[1] <=3) { print $r[0]."\n"; }' < $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.frames > $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.dont_flipus

rm -f $READFILE_BASENAME.real.cons.fasta.dont_flipus.seq
for entry in `cat $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.dont_flipus`; do $BIN_DIR/fasta_header_grep.pl $entry$ $REAL_CONS_BASENAME.clustalnames.fasta >> $READFILE_BASENAME.real.cons.fasta.dont_flipus.seq ; done

cat $READFILE_BASENAME.real.cons.fasta.flipped $READFILE_BASENAME.real.cons.fasta.dont_flipus.seq > $READFILE_BASENAME.real.cons.fasta.flipped.seq

# merge the flipped phrap-arsfadig+ and the flipped real-arsfadig+

# begin new version merge

# in
cut -d'_' -f1 $READFILE_BASENAME.clustalnames.nonsingletons.pep.better_junk_in_phrap.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.better_junk_in_phrap.list.f1
$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.better_junk_in_phrap.list.f1 $READFILE_BASENAME.contigs.ARSFADIG.3mism.flipped.seq > $READFILE_BASENAME.clustalnames.nonsingletons.better_junk_in_phrap.flipped.seq

cut -d'_' -f1 $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.better_junk_after_realign.list > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.better_junk_after_realign.list.f1
$BIN_DIR/fasta_header_grep.pl -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.better_junk_after_realign.list.f1 $READFILE_BASENAME.real.cons.fasta.flipped.seq > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.better_junk_after_realign.flipped.seq

# out
cat $READFILE_BASENAME.clustalnames.nonsingletons.better_junk_in_phrap.flipped.seq $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.better_junk_after_realign.flipped.seq > $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq

# remove dbl2b contigs (RYSFADI, KYSFADI, RYSFADL)

grep \> $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq |sed -e 's/>//' >  $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.list

cut -f1 -d_ $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.dbl2b.list > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.dbl2b.list.f1
grep -w -v -f $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.stops_artificially_removed.dbl2b.list.f1 $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.list > $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.list.without_dbl2b

cp $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.with_dbl2b
$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.list.without_dbl2b  $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.with_dbl2b > $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq


# end new version merge

# get only phrap arsfadig+ with stops from the real bunch
#$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.1 $READFILE_BASENAME.real.cons.fasta.flipped.seq > $READFILE_BASENAME.real.cons.fasta.flipped.wsipc.seq

#cut -d'_' -f1 $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list.f1
#$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list.f1 $READFILE_BASENAME.contigs.ARSFADIG.3mism.flipped.seq > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.flipped.seq

# get any realigned arsfadig+ not already in the "final flipped set"
#cut -d'_' -f1 $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.not_in_phrap_ok.list > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.not_in_phrap_ok.list.f1
#$BIN_DIR/fasta_header_grep.pl -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.not_in_phrap_ok.list.f1 $READFILE_BASENAME.real.cons.fasta.flipped.seq > $READFILE_BASENAME.real.cons.fasta.flipped.not_in_phrap_ok.seq

#cat $READFILE_BASENAME.real.cons.fasta.flipped.not_in_phrap_ok.seq  $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.flipped.seq > $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq

# blast-search against our own dbl1-db

;;
( "blast" )
#cd $DATADIR; formatdb -p F -i pf_gbk_plasmodb_040512_var_gene_sequences.fasta ;cd -
#blastall -i $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq -d $DATA_DIR/pf_gbk_plasmodb_040512_var_gene_sequences.fasta -F F -o $READFILE_BASENAME"_vs_pf_gbk_040512.bln" -p blastn

BLN_FILE=$READFILE_BASENAME"_vs_"$BLASTDB".bln"
blastall -i $READFILE_BASENAME.clustalnames.contigs -d $DATA_DIR/$BLASTDB -F F -o $BLN_FILE -p blastn

$BLASTGFF=$BLN_FILE".c40i96.gff"
$BIN_DIR/parselargemegablast.dev.pl -c 0.4 -i 0.96 < $BLN_FILE > $BLASTGFF

# count fractions
$BIN_DIR/print_isolate_list_from_ace_file.pl --isolatefrac < $ACEFILE > $ACEFILE.isolatefrac
$BIN_DIR/print_isolate_list_from_ace_file.pl < $ACEFILE > $ACEFILE.frac_count
sort -k3,3 -k4,4nr $ACEFILE.frac_count > $ACEFILE.frac_count.sort

;;
( "tree" )
#
# msa & tree construction
#

# clustalw nt tree
clustalw $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq

# clustalw pept tree
perl -ne 'chomp; if (m/^>.*(Contig\d+)/) { print ">",$1,"\n"; } else { print $_,"\n"; }' < $READFILE_BASENAME.par.clustalnames.nonsingletons.pep > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.clustalnames
#

if [ $3 == "blosum" ] ;
then 
    echo "Using a BLOSUM matrix for clustal alignment for full tree."
    clustalw -INFILE=$READFILE_BASENAME.par.clustalnames.nonsingletons.pep.clustalnames -MATRIX=BLOSUM -PWMATRIX=BLOSUM 
else 
    echo "Using a GONNET matrix for clustal alignment for full tree."
    clustalw -INFILE=$READFILE_BASENAME.par.clustalnames.nonsingletons.pep.clustalnames -MATRIX=BLOSUM -PWMATRIX=BLOSUM 
fi


# count dominance and get most dominant
#$BIN_DIR/print_useful_info_about_each_contig.pl --typesort $ACEFILE.isolatefrac --fractioncount $ACEFILE.frac_count.sort --dominant_cutoff 5 --blasthitgff $READFILE_BASENAME"_vs_pf_gbk_040512.bln.c40i96.gff" --exclude_from_dominance $DOMINANCE_EXCLUDE --dominants_only > $READFILE_BASENAME.5dominants.info
$BIN_DIR/print_useful_info_about_each_contig.pl --typesort $ACEFILE.isolatefrac --fractioncount $ACEFILE.frac_count.sort --dominant_cutoff $DOMINANT_CUTOFF --blasthitgff $BLASTGFF --exclude_from_dominance $DOMINANCE_EXCLUDE --dominants_only > $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info
cut -f1  $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info > $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info.list
$BIN_DIR/fasta_header_grep.pl -w -f $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info.list $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.clustalnames > $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep

# add useful info to the trees

$BIN_DIR/add_useful_info_to_tree.pl --typesort $ACEFILE.isolatefrac --fractioncount $ACEFILE.frac_count.sort --dominant_cutoff $DOMINANT_CUTOFF --blasthitgff $BLASTGFF --tree $READFILE_BASENAME.clustalnames.nonsingletons.flipped.dnd --exclude_from_dominance $DOMINANCE_EXCLUDE > $READFILE_BASENAME.clustalnames.nonsingletons.flipped.info.dnd

$BIN_DIR/add_useful_info_to_tree.pl --typesort $ACEFILE.isolatefrac --fractioncount $ACEFILE.frac_count.sort --dominant_cutoff $DOMINANT_CUTOFF --blasthitgff $BLASTGFF --tree $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.dnd --exclude_from_dominance $DOMINANCE_EXCLUDE > $READFILE_BASENAME.par.clustalnames.nonsingletons.info.dnd

#
# branch merges 
#

#echo -ne "\nChecking for too narrow branch clusters in nt tree\n" >> $STATS

$BIN_DIR/show_clusters.pl $READFILE_BASENAME.clustalnames.nonsingletons.flipped.dnd 0.02 > $READFILE_BASENAME.clustalnames.nonsingletons.flipped.dnd.bushes 
#grep -c New $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.bushes >> $STATS
#$BIN_DIR/show_clusters.pl 0.02 <  $READFILE_BASENAME.par.clustalnames.nonsingletons.5dominants.pep  
#$BIN_DIR/show_clusters.pl 0.02 <  $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.dnd

;;
( "motifs" )
#
# find motifs
#

# fairly slow for full alignment -- use 5-most-dominants instead!
if ! [ -e $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep ] ;
then 
    perl -ne 'chomp; if (m/^>.*(Contig\d+)/) { print ">",$1,"\n"; } else { print $_,"\n"; }' < $READFILE_BASENAME.par.clustalnames.nonsingletons.pep > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.clustalnames
    $BIN_DIR/print_useful_info_about_each_contig.pl --typesort $ACEFILE.isolatefrac --fractioncount $ACEFILE.frac_count.sort --dominant_cutoff $DOMINANT_CUTOFF --blasthitgff $BLASTGFF --exclude_from_dominance $DOMINANCE_EXCLUDE --dominants_only > $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info
    cut -f1  $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info > $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info.list
    $BIN_DIR/fasta_header_grep.pl -w -f $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info.list $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.clustalnames > $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep
fi

if [ $3 == "blosum" ] ; 
then 
    echo "Using a BLOSUM matrix for clustal alignment for motif finding."
    clustalw -INFILE=$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep -MATRIX=BLOSUM -PWMATRIX=BLOSUM
else 
    echo "Using a GONNET matrix for clustal alignment for motif finding."
    clustalw -INFILE=$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep
fi

#$BIN_DIR/add_useful_info_to_tree.pl --typesort $ACEFILE.isolatefrac --fractioncount $ACEFILE.frac_count.sort --dominant_cutoff 6 --blasthitgff $BLASTGFF --tree $READFILE_BASENAME.par.clustalnames.nonsingletons.5dominants.dnd --exclude_from_dominance $DOMINANCE_EXCLUDE > $READFILE_BASENAME.par.clustalnames.nonsingletons.5dominants.info.dnd

mkdir $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.motifs
cd $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.motifs
# local copy of aln file to simplify output redirection
mkdir tmp
cp ../$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.aln .
#$BIN_DIR/section_alignent.pl --debug --in $READFILE_BASENAME.par.clustalnames.nonsingletons.5dominants.aln --fractioncount ../$ACEFILE.frac_count.sort --exclude_from_dominance $DOMINANCE_EXCLUDE --blasthitgff ../$BLASTGFF --dominant_cutoff 5 > ../$READFILE_BASENAME.par.clustalnames.nonsingletons.pep.5dominants.aln.section 2>&1
TEMPDIR=`pwd`/tmp && $BIN_DIR/section_alignent.pl --debug --in $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.aln --fractioncount ../$ACEFILE.frac_count.sort --exclude_from_dominance $DOMINANCE_EXCLUDE --blasthitgff ../$BLASTGFF --dominant_cutoff $DOMINANT_CUTOFF --countmodel oneperisolate > ../$READFILE_BASENAME.par.clustalnames.nonsingletons.pep.$DOMINANT_CUTOFFdominants.aln.opi.section 2>&1
# remove local copy
# rm $READFILE_BASENAME.par.clustalnames.nonsingletons.5dominants.aln

# make them searchable
echo Make them searchable..

DB_NAME=$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.aln.db.fasta
makemat -P $DB_NAME -d $READFILE_BASENAME
copymat -P $DB_NAME -r F
cd ..

exit

# Note - below is only for large set!

echo Section alignment..
mkdir $READFILE_BASENAME.par.clustalnames.nonsingletons.motifs
cd $READFILE_BASENAME.par.clustalnames.nonsingletons.motifs
cp ../$READFILE_BASENAME.par.clustalnames.nonsingletons.pep.aln .
# local copy of aln file to simplify output redirection

mkdir tmp
TEMPDIR=`pwd`/tmp && $BIN_DIR/section_alignent.pl --debug --in $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.aln --fractioncount ../$ACEFILE.frac_count.sort --exclude_from_dominance $DOMINANCE_EXCLUDE --blasthitgff ../$BLASTGFF --dominant_cutoff 7 --countmodel oneperisolate > ../$READFILE_BASENAME.par.clustalnames.nonsingletons.pep.aln.opi.section 2>&1

echo Make motifs searchable..
DB_NAME=$READFILE_BASENAME.par.clustalnames.nonsingletons.pep.aln.db.fasta
makemat -P $DB_NAME -d $READFILE_BASENAME
copymat -P $DB_NAME -r F
cd ..

exit

# calculate prob of drawing the combos (?)
# validate/provide statistics

# input test set as query
# parse results, esp connect motif P-vals & hit E-vals
echo Scoretest..
cd scoretest
for file in *fasta ; do 
    # -h 0.005 default for impala eval multipass inclusion threshold
    impala -P ../$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.motifs/$DB_NAME -i $file -o $file.imp -F F
    $BIN_DIR/score_new_sequences.pl --imp $file.imp --section ../$READFILE_BASENAME.par.clustalnames.nonsingletons.pep.$DOMINANT_CUTOFFdominants.aln.opi.section > $file.imp.scores
done
cd ..

;;
( "lmo" )
#
# split set into test and training sets
#

DB_NAME=$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.aln.db.fasta

cd scoretest

mkdir subdiv1 subdiv2 subdiv3 subdiv4 subdiv5
for dir in subdiv* ; do 
    cd $dir
    cp ../../$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.motifs/$DB_NAME .
    mv $DB_NAME $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep
    $BIN_DIR/split_fasta_into_sets.pl $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep

#        clustalw $READFILE_BASENAME.par.clustalnames.nonsingletons.5dominants.pep.p.set
    mv $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep.p.set p.set.$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep
    mv $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep.q.set q.set.$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep
    ln -s ../../$ACEFILE.frac_count.sort p.set.$ACEFILE.frac_count.sort
    ln -s ../../$ACEFILE.isolatefrac p.set.$ACEFILE.isolatefrac
    ln -s ../../$BLASTGFF p.set.$BLASTGFF
        
    $BIN_DIR/process_reads.sh p.set.$READFILE_BASENAME motifs

#    rm $DB_NAME
    mkdir p
    mkdir q
    cd p
    seqretsplit -sequence ../p.set.$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep -outseq dummy.fasta
    
    cd ../q
    seqretsplit -sequence ../q.set.$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep  -outseq dummier.fasta

    for file in *fasta ; do impala -P ../p.set.$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.motifs/p.set.$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.aln.db.fasta -i $file -o $file.imp -F F; $BIN_DIR/score_new_sequences.pl --imp $file.imp --section ../p.set.$READFILE_BASENAME.par.clustalnames.nonsingletons.pep.$DOMINANT_CUTOFFdominants.aln.opi.section > $file.imp.scores; done
    
    grep overrep *scores |sort -k1.7n > score
    perl -ne 's/contig/Contig/; print' < score > score.nfix
    $BIN_DIR/add_useful_info_to_tree.pl --typesort ../../../$ACEFILE.isolatefrac --fractioncount ../../../$ACEFILE.frac_count.sort --dominant_cutoff $DOMINANT_CUTOFF --blasthitgff ../../../$BLASTGFF --exclude_from_dominance $DOMINANCE_EXCLUDE --tree score.nfix > score.info
    cd ..
    cd ..   
done


# score the testsets
#for dir in subdiv* ; do 
    
# construct motifs from 
#    clustalw 
#    $BIN_DIR/process_reads.sh motifs

#    for file in *fasta ; do
#	impala -P ../../../$READFILE_BASENAME.par.clustalnames.nonsingletons.5dominants.motifs/$DB_NAME -i $file -o $file.imp -F F
#	$BIN_DIR/score_new_sequences.pl --imp $file.imp --section ../../../$READFILE_BASENAME.par.clustalnames.nonsingletons.pep.5dominants.aln.opi.section > $file.imp.scores
#    done

#done


#
# initial questions: do most severe dominant contigs have at least one "severe" domain according to current scoring?
# do most severe contigs score as more severe than most other contifs? how many?




# write manuscript
# submit
exit

;;
( "patientlmo" )

$BIN_DIR/split_full_seq_into_test_sets.pl -s $READFILE_BASENAME".screen.qp.plq.seq" -n 5
$BIN_DIR/make_new_impala_matrices_for_scoring_scheme.pl

for file in *fasta ; do impala -P ../p.set.$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.motifs/p.set.$READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.aln.db.fasta -i $file -o $file.imp -F F; $BIN_DIR/score_new_sequences.pl --imp $file.imp --section ../p.set.$READFILE_BASENAME.par.clustalnames.nonsingletons.pep.$DOMINANT_CUTOFFdominants.aln.opi.section > $file.imp.scores; done

$BIN_DIR/score_new_sequences_w_subregion_filter.pl
$BIN_DIR/query_min_per_patient_with_order.pl

exit
;;
( "stats" )

# STATS=/tmp/tmp.log
# $READFILE_BASENAME.seq
# echo -ne "\n" 

echo -ne "\nNumber of retained high quality reads: " >> $STATS
grep -c \> $READFILE_BASENAME".screen.qp.plq.seq" >> $STATS

echo -ne "\nAssembly criteria: " >> $STATS
grep -A3 WA $ACEFILE |tail +2 >> $STATS

echo -ne "\nNumber of singlets: " >> $STATS
grep -c \> $ACEBASE".singlets" >> $STATS

echo -ne "Number of effective singlets in assembly - one read contigs: " >> $STATS
grep Contig $ACEFILE | cut -f2,4 -d' ' | sort -k2n | perl -ne 'if(m/\ 1$/) { print } else { }' |wc -l >> $STATS
echo -ne "Conversely, the number of non-singleton contigs: " >> $STATS
wc -l $READFILE_BASENAME.nonsingletons.list| $BIN_DIR/wscut.pl -f 1 >> $STATS
echo -ne "Reads in larger contigs: " >> $STATS
grep Contig $ACEFILE | cut -f2,4 -d' ' | awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

echo -ne "\nNumber of $ARSFADIG_MOTIF containing, nonsingleton contigs at 3 mismatches: " >> $STATS
wc -l $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list | $BIN_DIR/wscut.pl -f 1 >> $STATS
echo -ne "Counting translated peptides: " >> $STATS 
grep -c \> $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep >> $STATS
echo -ne "Number of reads in such $ARSFADIG_MOTIF containing contigs: " >> $STATS
cut -d_ -f1 $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.frameless
grep -w -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.frameless $ACEFILE|cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

echo -ne "\nNumber of $ARSFADIG_MOTIF-negative contigs: " >> $STATS
wc -l $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.hitless |$BIN_DIR/wscut.pl -f 1 >> $STATS
echo -ne "Number of $ARSFADIG_MOTIF-negative reads: " >> $STATS
grep -w -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.frameless $ACEFILE|cut -f2,4 -d' '| awk '($2 > 1) {print}'|cut -f1 -d' '|sed -e's/Contig//'|perl -e '$current_id=1; while($c = <STDIN>) { chomp $c; while ( $current_id < $c) { print "Contig$current_id\n"; $current_id++ } $current_id++}' > $READFILE_BASENAME.clustalnames.ARSFADIG.3mism.list.hitless
grep -w -f $READFILE_BASENAME.nonsingletons.list $READFILE_BASENAME.clustalnames.ARSFADIG.3mism.list.hitless > $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.hitless
grep -w -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.hitless $ACEFILE|cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

# count stop codons and find out how many of the MOTIFless contigs have
echo -ne "Number of $ARSFADIG_MOTIF-less contigs with at least one stop codon free orf: " >> $STATS
$BIN_DIR/fasta_header_grep.pl -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.hitless $ACEBASE.contigs >  $READFILE_BASENAME.nonsingletons.ARSFADIG.3mism.list.hitless.fasta
transeq -frame 6 -sequence $READFILE_BASENAME.nonsingletons.ARSFADIG.3mism.list.hitless.fasta -outseq $READFILE_BASENAME.nonsingletons.ARSFADIG.3mism.list.hitless.pep

perl -ne 'chomp; if (m/^>.+(Contig\d+_\d+)/) { print ">",$1,"\n"; } else { print $_,"\n"; }' < $READFILE_BASENAME.nonsingletons.ARSFADIG.3mism.list.hitless.pep > $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless.pep

$BIN_DIR/fasta_pep_re.pl --regexp \\* --infile $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless.pep |cut -f3 |grep Contig |sort -k1.7n |uniq -c |sort -n > $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless.pep.nstopsort

$BIN_DIR/wscut.pl -f 3 < $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless.pep.nstopsort |sort > $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless.pep.nstopsort.list
cut -f1 -d'_' $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless.pep.nstopsort.list |uniq > $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless.pep.nstopsort.frameless

grep -w -v -f $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless.pep.nstopsort.frameless $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless > $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.nostop.list
wc -l $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.nostop.list |$BIN_DIR/wscut.pl -f 1 >> $STATS

# with 1 stop?
echo -ne "Number of $ARSFADIG_MOTIF-less contigs with at least one ORF with at most one stop: " >> $STATS
grep \ 1\  $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless.pep.nstopsort |perl -ne 'm/(Contig\d+)_/; print $1,"\n";' |sort|uniq > $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless.pep.onestop
wc -l $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless.pep.onestop |$BIN_DIR/wscut.pl -f 1 >> $STATS

echo -ne "Number of reads: " >> $STATS
grep -w -f $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.list.hitless.pep.onestop $ACEFILE|cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

echo -ne "\nNumber of stop codon containing reading frames from $ARSFADIG_MOTIF positive contigs: " >> $STATS
wc -l $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort |$BIN_DIR/wscut.pl -f 1 >> $STATS

echo -ne "Nr of reads in these stop codon containing $ARSFADIG_MOTIF positive contigs: "  >> $STATS
cut -d_ -f1 $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.frameless
grep -w -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.frameless $ACEFILE|cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

echo -ne "Number of stop codon free $ARSFADIG_MOTIF positive reading frames: " >> $STATS
grep -c \> $READFILE_BASENAME.clustalnames.nonsingletons.nostop.pep >> $STATS

echo -ne "Number of reads in stop codon free $ARSFADIG_MOTIF positive reading frames: " >> $STATS
cut -d_ -f1 $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list > $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list.frameless
grep -w -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list.frameless $ACEFILE|cut -f2,4 -d' '|awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

echo -ne "\nNumber of RE-AL stop codon containing $ARSFADIG_MOTIF positive reading frames: " >> $STATS
wc -l $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort |$BIN_DIR/wscut.pl -f 1 >> $STATS

echo -ne "Number of reads in RE-AL stop codon containing $ARSFADIG_MOTIF positive reading frames: " >> $STATS
cut -d_ -f1 $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.frameless
grep -w -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.frameless $ACEFILE |cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

echo -ne "Number of RE-AL stop codon free $ARSFADIG_MOTIF positive reading frames: " >> $STATS
grep -c \> $REAL_CONS_BASENAME.clustalnames.nonsingletons.nostop.pep >> $STATS

echo -ne "Number of reads in RE-AL stop codon free $ARSFADIG_MOTIF positive contigs: " >> $STATS
cut -d_ -f1 $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nostop.list > $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nostop.list.frameless
grep -w -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nostop.list.frameless $ACEFILE |cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

# nr of 1 contigs merged from phrap
echo -ne "\nNumber of original phrap contigs retained in merge: " >> $STATS
wc -l $READFILE_BASENAME.clustalnames.nonsingletons.pep.stopcount.better_junk_in_phrap |$BIN_DIR/wscut.pl -f 1 >> $STATS

echo -ne "Number of reads: " >> $STATS
grep -w -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.stopcount.better_junk_in_phrap $ACEFILE |cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

# merged from real
echo -ne "Number of re-aligned phrap contigs retained in merge: " >> $STATS
wc -l $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.stopcount.better_junk_after_realign |$BIN_DIR/wscut.pl -f 1 >> $STATS

echo -ne "Number of reads: " >> $STATS
grep -w -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.stopcount.better_junk_after_realign $ACEFILE |cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

# nr of 1-stop contigs retained at this stage?

echo -ne "\nNumber of FINAL SET (1) stop codon containing $ARSFADIG_MOTIF positive reading frames: " >> $STATS
wc -l $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort |$BIN_DIR/wscut.pl -f 1 >> $STATS
echo -ne "Number of reads in those contigs: " >> $STATS

$BIN_DIR/wscut.pl -f 3 < $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort |sort > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort.list
cut -d_ -f1 $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort.list > $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort.list.frameless
grep -w -f $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort.list.frameless $ACEFILE|cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

echo -ne "\nNumber of top $DOMINANT_CUTOFF-dominant contigs: " >> $STATS
wc -l $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info.list |$BIN_DIR/wscut.pl -f 1 >> $STATS
echo -ne "Number of reads: " >> $STATS
grep -w -f $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info.list $ACEFILE|cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

echo -ne "\nNumber of top 3-af-top-dominant contigs: " >> $STATS
wc -l $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info.list |$BIN_DIR/wscut.pl -f 1 >> $STATS
echo -ne "Number of reads: " >> $STATS
grep -w -f $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info.list $ACEFILE|cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

echo -ne "\nNumber of top $DOMINANT_CUTOFF-dominant contigs in final contig set: " >> $STATS
grep -c \> $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep >> $STATS
echo -ne "\nNumber of included reads: " >> $STATS
grep \> $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.pep |sed -e 's/>//'> $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.list
grep -w -f $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.list $ACEFILE|cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

echo -ne "Number of top $DOMINANT_CUTOFF-dominant contigs in non-$ARSFADIG_MOTIF contigs: " >> $STATS
grep -w -f $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info.list $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.hitless >  $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.hitless.$DOMINANT_CUTOFFdominants
wc -l $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.hitless.$DOMINANT_CUTOFFdominants |$BIN_DIR/wscut.pl -f 1 >> $STATS
echo -ne "Number of reads: " >> $STATS
grep -w -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.hitless.$DOMINANT_CUTOFFdominants $ACEFILE|cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}' >> $STATS

echo -ne "Number of top $DOMINANT_CUTOFF-dominant contigs neither in final set nor in non-$ARSFADIG_MOTIF contifs: " >> $STATS
grep -v -f $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.list $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info.list > $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.not_in_final
grep -v -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.hitless.$DOMINANT_CUTOFFdominants $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.not_in_final > $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.not_in_final.minus_hitless
wc -l $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.not_in_final.minus_hitless |$BIN_DIR/wscut.pl -f 1 >> $STATS
echo -ne "Number of reads: " >> $STATS
grep -w -f $READFILE_BASENAME.par.clustalnames.nonsingletons.$DOMINANT_CUTOFFdominants.not_in_final.minus_hitless $ACEFILE|cut -f2,4 -d' '| awk 'BEGIN {count=0;} ($2 > 1) {count=count+$2;} END {print count}'  >> $STATS

echo -ne "\nChecking for too narrow branch clusters in nt tree\n" >> $STATS

#$BIN_DIR/show_clusters.pl $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq 0.02 > $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.bushes 
#grep -c New $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.bushes >> $STATS


#if [ -e $READFILE_BASENAME.par.clustalnames.nonsingletons.5dominants.pep ] ;
#then 
#fi

;;
( "cleanup" )
#
# clean up
#

rm -f $READFILE_BASENAME.nonsingletons.list $READFILE_BASENAME.nonsingletons.fasta $READFILE_BASENAME.nonsingletons.pep $READFILE_BASENAME.nonsingletons.clustalnames.pep
rm -f $READFILE_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed $READFILE_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.fuzzpro 

rm -f $READFILE_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.list.unedited $READFILE_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes.1 
rm -f $READFILE_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes
rm -f $READFILE_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes.list $READFILE_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list 
rm -f $READFILE_BASENAME.nonsingletons.clustalnames.ARSFADIG.3mism.pep $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep
rm -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list
rm -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list $READFILE_BASENAME.clustalnames.nonsingletons.nostop.pep $READFILE_BASENAME.clustalnames.nonsingletons.nostop.pep
rm -f $READFILE_BASENAME.real.cons.fasta $REAL_CONS_BASENAME.clustalnames.fasta $REAL_CONS_BASENAME.nonsingletons.clustalnames.fasta 
rm -f $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.fuzzpro 
rm -f $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.stops_artificially_removed.ARSFADIG.3mism.list.unedited $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes.1
rm -f $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes $REAL_CONS_BASENAME.nonsingletons.clustalnames.pep.ARSFADIG.3mism.list.unedited.dupes.list
rm -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.pep $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort

rm -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nostop.list $REAL_CONS_BASENAME.clustalnames.nonsingletons.nostop.pep $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.1
rm -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.with_stops_in_phrap_cons $READFILE_BASENAME.par.clustalnames.nonsingletons.pep 
rm -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list

rm -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.not_in_phrap_ok.pep $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.not_in_phrap_ok.list 
rm -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.list.f1 $REAL_CONS_BASENAME.fasta.flipus.seq $REAL_CONS_BASENAME.fasta.flipped.seq $REAL_CONS_BASENAME.fasta.flipped
rm -f $REAL_CONS_BASENAME.fasta.dont_flipus.seq  $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.dont_flipus $REAL_CONS_BASENAME.fasta.flipped.not_in_phrap_ok.seq 
rm -f $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.ARSFADIG.3mism.not_in_phrap_ok.list.f1
rm -f $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.frames $READFILE_BASENAME.clustalnames.contigs $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.flipus $READFILE_BASENAME.contigs.ARSFADIG.3mism.clustalnames.flipus.seq $READFILE_BASENAME.contigs.ARSFADIG.3mism.flipped $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.dont_flipus $READFILE_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.dont_flipus.seq $READFILE_BASENAME.contigs.ARSFADIG.3mism.flipped.seq $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.frames $REAL_CONS_BASENAME.clustalnames.nonsingletons.ARSFADIG.3mism.flipus
rm -f $READFILE_BASENAME.par.clustalnames.nonsingletons.pep $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq

rm -f $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.aln $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.aln
rm -f $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.dnd $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.dnd

rm -f $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.aln.section.debug $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.info.dnd $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.info.dnd 
rm -f $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info
rm -rf $READFILE_BASENAME.raws $READFILE_BASENAME.real $READFILE_BASENAME.real.fasta $READFILE_BASENAME.aln.real $READFILE_BASENAME.cons

rm -f $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort.list.frameless  $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nostop.list.frameless
rm -f  $REAL_CONS_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.frameless  $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.list.frameless  $READFILE_BASENAME.clustalnames.nonsingletons.pep.nstopsort.list.frameless

rm -f $READFILE_BASENAME.clustalnames.nonsingletons.pep.nostop.flipped.seq $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.nstopsort.list
rm -f $READFILE_BASENAME.clustalnames.nonsingletons.flipped.seq.bushes $READFILE_BASENAME.par.clustalnames.nonsingletons.info.dnd $READFILE_BASENAME.$DOMINANT_CUTOFFdominants.info.list 
rm -f $READFILE_BASENAME.clustalnames.nonsingletons.flipped.aln $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.clustalnames $READFILE_BASENAME.clustalnames.nonsingletons.flipped.dnd
rm -f $READFILE_BASENAME.clustalnames.nonsingletons.flipped.info.dnd $ACEFILE.frac_count $ACEFILE.frac_count.sort $ACEFILE.isolatefrac 

rm -f  $BLN_FILE $BLASTGFF

cp $STATS $STATS.old 
rm -f $STATS

rm -f motifs $READFILE_BASENAME.par.clustalnames.nonsingletons.pep.$DOMINANT_CUTOFFdominants.aln.section

;;
( "most" ) 
$BIN_DIR/process_reads.sh $READFILE translate
$BIN_DIR/process_reads.sh $READFILE realign
$BIN_DIR/process_reads.sh $READFILE flip
$BIN_DIR/process_reads.sh $READFILE blast
$BIN_DIR/process_reads.sh $READFILE tree
$BIN_DIR/process_reads.sh $READFILE motifs
$BIN_DIR/process_reads.sh $READFILE stats


;;
( "all" )
$BIN_DIR/process_reads.sh $READFILE trim
$BIN_DIR/process_reads.sh $READFILE phrap
$BIN_DIR/process_reads.sh $READFILE translate
$BIN_DIR/process_reads.sh $READFILE realign
$BIN_DIR/process_reads.sh $READFILE flip
$BIN_DIR/process_reads.sh $READFILE blast
$BIN_DIR/process_reads.sh $READFILE tree
$BIN_DIR/process_reads.sh $READFILE motifs
$BIN_DIR/process_reads.sh $READFILE stats
$BIN_DIR/process_reads.sh $READFILE cleanup

;;
esac
