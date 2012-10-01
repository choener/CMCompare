#!/bin/bash

# compare our results with infernal!

tail -n $1 data/output/nosamesame | while read line ; do
  #extract info from dataset
  SCORE=`echo "$line" | awk '{print $1}'`
  CM1=`echo "$line" | awk '{print $2}'`
  CM2=`echo "$line" | awk '{print $3}'`
  CM1=`basename $CM1`
  CM2=`basename $CM2`

  #redo calculation with rnastring output
  FULL=`./hsCMCompare Rfam/$CM1 Rfam/$CM2`
  STR=`echo $FULL | awk '{print $5}'`
  NEWSC=`echo $FULL | awk '{print $3 " " $4}'`

  #put string into temp file
  echo "> $CM1-$CM2" > data/infernal.tmp
  echo "$STR" >> data/infernal.tmp

  #run infernal
  IN1=`~/snorna-tools/bin/cmsearch --no-null3 --cyk --fil-no-hmm --fil-no-qdb Rfam/$CM1 data/infernal.tmp`
  IN2=`~/snorna-tools/bin/cmsearch --no-null3 --cyk --fil-no-hmm --fil-no-qdb Rfam/$CM2 data/infernal.tmp`

  #grep scores
  SIN1=`echo $IN1 | grep -o "Score = [0-9\.]*" | head -n1`
  SIN2=`echo $IN2 | grep -o "Score = [0-9\.]*" | head -n1`

  #put everything into an output file
  OUT=`basename $CM1 .cm`-`basename $CM2 .cm`.out
  echo "CMCOMPARE $SCORE $CM1 $CM2      $NEWSC" > data/infernal/$OUT
  echo "INFERNAL-1 $SIN1" >> data/infernal/$OUT
  echo "INFERNAL-2 $SIN2" >> data/infernal/$OUT
  echo "$STR" >> data/infernal/$OUT
  echo "" >> data/infernal/$OUT
  echo "" >> data/infernal/$OUT
  echo "" >> data/infernal/$OUT
  #and the whole output, too
  echo "$IN1" >> data/infernal/$OUT
  echo "$IN2" >> data/infernal/$OUT
done
