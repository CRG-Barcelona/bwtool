#!/bin/bash

name=`basename $0 .sh`
./core-test.sh $name \
  answers/${name}.wig \
  tested.bw \
  1 fix no \
  wigs/main.wig \
  ../../bwtool remove less 3 -decimals=1 main.bw tested.bw
exit $?
