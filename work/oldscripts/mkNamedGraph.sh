#! /bin/bash

repWithName () {
  NUM1=`echo $2 | grep -o "[0-9]*"`
  NUM2=`echo $3 | grep -o "[0-9]*"`
  NAME1=`cat Rfam/RF$NUM1.cm | grep NAME | awk '{print $2}'`
  NAME2=`cat Rfam/RF$NUM2.cm | grep NAME | awk '{print $2}'`
  W=`echo $1 | awk '{printf "%d", $1}'`
  echo "\"$NUM1\n$NAME1\" -- \"$NUM2\n$NAME2\" [label = \"$W\"];"
}

echo "graph g {"

tail -n$1 data/output/nosamesame | while read line ; do
  repWithName $line
done | sort
# | awk '{printf "\"%s\" -- \"%s\" [label = \"%d\"];\n",substr($2,12,5),substr($3,12,5),$1}' | sort

echo "}"
