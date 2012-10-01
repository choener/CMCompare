echo "graph g {"
tail -n$1 data/output/nosamesame | awk '{printf "\"%s\" -- \"%s\" [label = \"%d\"];\n",substr($2,12,5),substr($3,12,5),$1}' | sort
echo "}"
