#! /bin/bash

cms=( $(ls -Sr testdata/*cm) )
len=$((${#cms[@]} -1))

for i in $(seq 0 $len)
  do
    for j in $(seq $i $len)
    do
      echo "${cms[$i]}" "${cms[$j]}"
    done
  done
