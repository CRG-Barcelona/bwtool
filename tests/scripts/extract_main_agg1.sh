#!/bin/bash

name=`basename $0 .sh`
./core-test.sh $name \
  answers/${name}.bed \
  tested.bed \
  0 0 0 \
  wigs/main.wig \
  ../../bwtool extract bed ../beds/agg1.bed main.bw tested.bed
exit $?
