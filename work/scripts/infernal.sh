#!/bin/bash

# compare our results with infernal!
# args:
# 1 : cm 1
# 2 : cm 2
# 3 : output dir

#extract info from dataset
CM1=`basename $1`
CM2=`basename $2`

OUT=`basename $CM1 .cm`-`basename $CM2 .cm`.out
if [ ! -f $3/$OUT ] ;
then
  #redo calculation with rnastring output
  FULL=`./CMCompare $1 $2`
  STR=`echo $FULL | awk '{print $5}'`
  SCORES=`echo $FULL | awk '{print $3 "   " $4}'`

  #put string into temp file
  echo "> $CM1-$CM2" > $3/tmp.tmp
  echo "$STR" >> $3/tmp.tmp

  #run infernal
  IN1=`~/bin/cmsearch --no-qdb --no-null3 --cyk --fil-no-hmm --fil-no-qdb $1 $3/tmp.tmp`
  IN2=`~/bin/cmsearch --no-qdb --no-null3 --cyk --fil-no-hmm --fil-no-qdb $2 $3/tmp.tmp`

  #grep scores
  SIN1=`echo $IN1 | egrep -o "Score = (-)?[0-9\.]*" | head -n1`
  SIN2=`echo $IN2 | egrep -o "Score = (-)?[0-9\.]*" | head -n1`

  #put everything into an output file
  echo "CMCOMPARE $CM1   $CM2 === $SCORES === $SIN1   $SIN2" > $3/$OUT
  echo "$FULL" >> $3/$OUT
  echo "$STR" >> $3/$OUT
  echo "" >> $3/$OUT
  echo "" >> $3/$OUT
  echo "" >> $3/$OUT
  #and the whole output, too
  echo "$IN1" >> $3/$OUT
  echo "$IN2" >> $3/$OUT
fi
