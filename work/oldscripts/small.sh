#!/bin/bash

# compare our results with infernal!

tail -n $1 data/smallmodels-results/res | while read line ; do
  #extract info from dataset
  CM1=`echo "$line" | awk '{print $1}'`
  CM2=`echo "$line" | awk '{print $2}'`
  CM1=`basename $CM1`
  CM2=`basename $CM2`

  #redo calculation with rnastring output
  FULL=`./hsCMCompare Rfam/$CM1 Rfam/$CM2`
  STR=`echo $FULL | awk '{print $5}'`
  NEWSC=`echo $FULL | awk '{print $3 " " $4}'`

  #put string into temp file
  echo "> $CM1-$CM2" > data/smallmodels-results.tmp
  echo "$STR" >> data/smallmodels-results.tmp

  #run infernal
  IN1=`~/bin/cmsearch --no-null3 --cyk --fil-no-hmm --fil-no-qdb Rfam/$CM1 data/smallmodels-results.tmp`
  IN2=`~/bin/cmsearch --no-null3 --cyk --fil-no-hmm --fil-no-qdb Rfam/$CM2 data/smallmodels-results.tmp`

  #grep scores
  SIN1=`echo $IN1 | egrep -o "Score = (-)?[0-9\.]*" | head -n1`
  SIN2=`echo $IN2 | egrep -o "Score = (-)?[0-9\.]*" | head -n1`

  #put everything into an output file
  OUT=`basename $CM1 .cm`-`basename $CM2 .cm`.out
  echo "CMCOMPARE $CM1 $CM2      $NEWSC  $SIN1   $SIN2" > data/smallmodels-results/$OUT
  echo "$STR" >> data/smallmodels-results/$OUT
  echo "" >> data/smallmodels-results/$OUT
  echo "" >> data/smallmodels-results/$OUT
  echo "" >> data/smallmodels-results/$OUT
  #and the whole output, too
  echo "$IN1" >> data/smallmodels-results/$OUT
  echo "$IN2" >> data/smallmodels-results/$OUT
done
