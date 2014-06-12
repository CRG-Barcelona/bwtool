#!/bin/bash

name=`basename $0 .sh`
./core-test.sh $name \
  answers/${name}.txt \
  tested.txt \
  0 0 0 \
  wigs/main.wig wigs/second.wig \
  ../../bwtool agg 3:3 ../beds/r0.bed,../beds/yellow.bed main.bw,second.bw tested.txt
exit $?
