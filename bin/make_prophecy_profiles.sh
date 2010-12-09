#!/bin/bash

SECTIONFILE=$1
DATABASE=$2
Efract=300

fasta_dump_for_matrices.pl --section $SECTIONFILE

# prepare for significant score estimates by making a scrambled database
shuffleseq -shuffle 100 $DATABASE $DATABASE.shuffled

for fastafile in  *fasta
do
  filebasename=`basename $fastafile .fasta`

  clustalw -INFILE=$fastafile -MATRIX=BLOSUM -PWMATRIX=BLOSUM -OUTFILE=$filebasename.aln

  prophecy -sequence $filebasename.aln -type G -name $filebasename -outfile $filebasename.prophecy -datafile EBLOSUM62 -open 3 -extension 0.3

  # significance check: scramble & score several times, note highest scores 

  prophet -sequence $DATABASE.shuffled -infile $filebasename.prophecy -outfile $DATABASE.$filebasename.shuffled.prophet -gapopen 3 -gapextend 0.3

  nDATABASEseq=`grep -c \> $DATABASE`
  nCutoff=$(($nDATABASEseq*100/$Efract));

  scoreCutoff=`sort_prophet_hits.pl --in $DATABASE.$filebasename.shuffled.prophet --topn $nCutoff |grep ^Score: |cut -d' ' -f2 |tail -1`
  echo "At $nCutoff (1/$Efract) expected random hits among $nDATABASEseq * 100 the score cutoff is $scoreCutoff."
  
  prophet -sequence $DATABASE -infile $filebasename.prophecy -outfile $DATABASE.$filebasename.prophet -gapopen 3 -gapextend 0.3
  sort_prophet_hits.pl --in $DATABASE.$filebasename.prophet --threshold $scoreCutoff >  $DATABASE.$filebasename.prophet.significant

  significantHits=`grep -c Local\: $DATABASE.$filebasename.prophet.significant`
  echo "Found $significantHits significant hits in $DATABASE with $filebasename at cutoff $scoreCutoff."
done

# order hits according to scores 
