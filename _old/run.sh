#! /bin/bash

cms=( $(ls testdata/*) )
len=$((${#cms[@]} -1))

for i in $(seq 0 $len)
  do
    for j in $(seq $i $len)
    do
      ./CMCompare --rna "${cms[$i]}" "${cms[$j]}"
    done
  done
