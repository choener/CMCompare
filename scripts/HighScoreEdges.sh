#! /bin/bash

# Christian Hoener zu Siederdissen, 2010
# choener@tbi.univie.ac.at
# Distributed "as-is" without any warranty

# Generate a graph of high-scoring edges between nodes. This script does not
# check parameters before using them. This is mostly because no operation
# writes or deletes anything. All results go to stdout!

# The following parameters are required:
# 1. directory containing Rfam models
# 2. number of edges to generate
# x. input is stdin, with (1) model name (2) model name (3) score

# Everything is then wrapped for graphviz dot output:
# zcat rfam-weakpairs.gz | HighScoreEdges.sh ../Rfam 200 | dot -Tps > graph.ps

repWithName () {
  NUM1=`echo $2 | grep -o "[0-9]*"`
  NUM2=`echo $3 | grep -o "[0-9]*"`
  NAME1=`cat $1/RF$NUM1.cm | grep NAME | awk '{print $2}'`
  NAME2=`cat $1/RF$NUM2.cm | grep NAME | awk '{print $2}'`
  W=`echo $4 | awk '{printf "%d", $1}'`
  echo "\"$NUM1\n$NAME1\" -- \"$NUM2\n$NAME2\" [label = \"$W\"];"
}

echo "graph g {"

tail -n$2 | while read line ; do
  repWithName $1 $line
done | sort

echo "}"
