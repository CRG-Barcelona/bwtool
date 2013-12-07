#!/bin/bash

name=`basename $0 .sh`
./core-test.sh $name \
  answers/${name}.txt \
  tested.txt \
  0 0 0 \
  wigs/main.wig \
  ../../bwtool agg -ends 3:3 ../beds/agg1.bed main.bw tested.txt
exit $?
