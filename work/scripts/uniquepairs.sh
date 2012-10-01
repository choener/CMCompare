#! /bin/bash

cms=( $(ls -Sr $1) )
len=$((${#cms[@]} -1))

for i in $(seq 0 $len)
  do
    for j in $(seq $(($i +1)) $len)
    do
      echo "$1/${cms[$i]}" "$1/${cms[$j]}"
    done
  done
