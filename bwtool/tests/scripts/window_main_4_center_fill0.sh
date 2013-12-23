#!/bin/bash

name=`basename $0 .sh`
./core-test.sh $name \
  answers/${name}.txt \
  tested.txt \
  0 0 0 \
  wigs/main.wig \
  ../../bwtool window 4 main.bw -o=tested.txt -center -fill=0
exit $?
