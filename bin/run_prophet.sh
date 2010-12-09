#!/bin/bash

#SECTIONFILE=$1
#fasta_dump_for_matrices.pl --section $SECTIONFILE

DATABASE=$2

fastafile=$3
filebasename=`basename $fastafile .fasta`

scoreCutoff=$4

# prepare for significant score estimates by making a scrambled database
#shuffleseq -shuffle 100 $DATABASE $DATABASE.shuffled

#clustalw -INFILE=$fastafile -MATRIX=BLOSUM -PWMATRIX=BLOSUM -OUTFILE=$filebasename.aln
#prophecy -sequence $filebasename.aln -type G -name $filebasename -outfile $filebasename.prophecy -datafile EBLOSUM62 -open 3 -extension 0.3

# significance check: scramble & score several times, note highest scores 
#prophet -sequence $DATABASE -infile $filebasename.prophecy -outfile $DATABASE.$filebasename.prophet -gapopen 3 -gapextend 0.3

#nDATABASEseq=`grep -c \> $DATABASE`

#nCutoff=$(($nDATABASEseq*100/$Efract));

#scoreCutoff=`sort_prophet_hits.pl --in $DATABASE.$filebasename.shuffled.prophet --topn $nCutoff |grep ^Score: |cut -d' ' -f2 |tail -1`
#echo "At $nCutoff (1/$Efract) expected random hits among $nDATABASEseq * 100 the score cutoff is $scoreCutoff."
  
prophet -sequence $DATABASE -infile $filebasename.prophecy -outfile $DATABASE.$filebasename.prophet -gapopen 3 -gapextend 0.3
sort_prophet_hits.pl --in $DATABASE.$filebasename.prophet --threshold $scoreCutoff >  $DATABASE.$filebasename.prophet.hits

significantHits=`grep -c Local\: $DATABASE.$filebasename.prophet.hits`
echo "Found $significantHits significant in $DATABASE with $filebasename over cutoff $scoreCutoff (see $DATABASE.$filebasename.prophet.hits)."
#cat $DATABASE.$filebasename.prophet.hits
# order hits according to scores 
