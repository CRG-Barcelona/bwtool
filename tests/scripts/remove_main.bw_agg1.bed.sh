#!/bin/bash

name=`basename $0 .sh`
./core-test.sh $name \
  answers/${name}.wig \
  tested.bw \
  0 var no \
  wigs/main.wig \
  ../../bwtool remove mask ../beds/agg1.bed -decimals=0 -wigtype=bg main.bw tested.bw
exit $?
