#!/bin/bash

name=`basename $0 .sh`
./core-test.sh $name \
  answers/${name}.txt \
  bad.txt \
  1 var no \
  wigs/main.wig \
  ../../bwtool lift main.bw ../misc/to_new_transposon_and_gap.chain lifted.bw -unlifted=bad.txt
exit $?
