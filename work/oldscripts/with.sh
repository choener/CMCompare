#!/bin/bash


exec 3< $1
while read <&3
do
  ./hsCMCompare $REPLY | tee -a ./output/`basename $1`.out
done
exec 3>&-
