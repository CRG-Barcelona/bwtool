#!/bin/bash

name=`basename $0 .sh`
./core-test.sh $name \
  answers/${name}.wig \
  lifted.bw \
  1 var no \
  wigs/main.wig \
  ../../bwtool lift main.bw ../misc/to_new_transposon_and_gap.chain lifted.bw
exit $?
