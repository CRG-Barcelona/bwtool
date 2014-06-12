#!/bin/bash

name=`basename $0 .sh`
./core-test.sh $name \
  answers/${name}.cg \
  tested.cg \
  0 0 0 \
  wigs/main.wig \
  ../../bwtool chromgraph -every=5 main.bw tested.cg
exit $?
