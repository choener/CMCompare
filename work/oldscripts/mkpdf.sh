#! /bin/bash
./mkNamedGraph.sh $1 | dot -Tpdf > data/graph.pdf
