#!/bin/bash

SECTIONFILE=$1
fasta_dump_for_matrices.pl --section $SECTIONFILE

DATABASE=$2
#Efract=300

# prepare for significant score estimates by making a scrambled database
if [ ! -e $DATABASE.shuffled ] ; 
then
    shuffleseq -shuffle 100 $DATABASE $DATABASE.shuffled ;
fi

fastafile=$3
filebasename=`basename $fastafile .fasta`

if [ ! -e $filebasename.aln ] ;
then
clustalw -INFILE=$fastafile -MATRIX=BLOSUM -PWMATRIX=BLOSUM -OUTFILE=$filebasename.aln
fi

prophecy -sequence $filebasename.aln -type G -name $filebasename -outfile $filebasename.prophecy -datafile EBLOSUM62 -open 3 -extension 0.3

  # significance check: scramble & score several times, note highest scores 

prophet -sequence $DATABASE.shuffled -infile $filebasename.prophecy -outfile $DATABASE.$filebasename.shuffled.prophet -gapopen 3 -gapextend 0.3

nDATABASEseq=`grep -c \> $DATABASE`

for Efract in {10,20,30,40,50,75,100,150,200,300,400,500,750,1000,10000}
do
  nCutoff=$(($nDATABASEseq*100/$Efract));
  
  scoreCutoff=`sort_prophet_hits.pl --in $DATABASE.$filebasename.shuffled.prophet --topn $nCutoff |grep ^Score: |cut -d' ' -f2 |tail -1`
  echo "At $nCutoff (1/$Efract) expected random hits among $nDATABASEseq * 100 the score cutoff is $scoreCutoff."
  
  prophet -sequence $DATABASE -infile $filebasename.prophecy -outfile $DATABASE.$filebasename.prophet -gapopen 3 -gapextend 0.3
  sort_prophet_hits.pl --in $DATABASE.$filebasename.prophet --threshold $scoreCutoff >  $DATABASE.$filebasename.prophet.significant
  
  significantHits=`grep -c Local\: $DATABASE.$filebasename.prophet.significant`
  echo "Found $significantHits significant hits in $DATABASE with $filebasename at cutoff $scoreCutoff."
done

# order hits according to scores 
