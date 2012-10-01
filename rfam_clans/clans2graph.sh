#!/bin/bash

# convert all clans to ingroup/outgroup graphs

ALLCLANS=`cat CLANDESC.all | egrep "^(AC|MB)" | awk '{print $2}'`

CNT=1

# prep clan files

for i in $ALLCLANS
do
  ARR[$CNT]=$i
  let CNT=$CNT+1
done

let CNT=$CNT-1

for i in `seq 1 $CNT`
do
  if [ `echo ${ARR[$i]} | head -c2` = "CL" ]
  then
    rm -f ${ARR[$i]}.tmp.clan
    j=$i+1
    while [ `echo ${ARR[$j]} | head -c2` = "RF" ]
    do
      RF=`echo ${ARR[$j]} | head -c7`
      echo $RF >> ${ARR[$i]}.tmp.clan
      let j=$j+1
      if [ $j -gt $CNT ]
      then
        break
      fi
    done
  fi
done

# from tmp.clan files, generate real files
for i in `ls *.tmp.clan`
do
  FNAME="`basename $i .tmp.clan`.clan"

  echo "graph g {" > $FNAME

  # now, we add all intra-clan edges
  for cm1 in `cat $i`
  do
    for cm2 in `cat $i`
    do
      if [[ $cm1 > $cm2 ]]
      then
        EDGE=`zcat rfam-weakpairs.gz | grep "$cm1" | grep "$cm2"`
        SCORE=`echo $EDGE | awk '{print $3}'`
        echo $SCORE
      fi
    done
  done

  echo "}" >> $FNAME
done
