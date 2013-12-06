#!/bin/bash

name=`basename $0 .sh`
./core-test.sh $name \
  answers/${name}.wig \
  tested.bw \
  1 var no \
  wigs/main.wig \
  ../../bwtool fill 0 main.bw tested.bw
exit $?
