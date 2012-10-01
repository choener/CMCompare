#! /bin/bash

cms=( $(ls input/input.*) )
len=$((${#cms[@]} -1))

#for i in `seq 0 5` ; do
#  for j in `seq 0 15` ; do
#    # running a job in the background in parallel:
#    # (ls &)
#    host="xc0$i"
#    job=(i*16+j)
#    echo $host
#    echo "${cms[$job]}"
#    (ssh $host "cd /scr/airline/choener/gogogo ; ./hsCMCompare" &)
#  done
#done

let i=$1

for j in `seq 0 7` ; do
  # running a job in the background in parallel:
  # (ls &)
  host="xc0$i"
  job=(i*8+j)
  echo $host
  echo "${cms[$job]}"
  (./with.sh "${cms[$job]}" &)
done
