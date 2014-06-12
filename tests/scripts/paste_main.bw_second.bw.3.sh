#!/bin/bash

name=`basename $0 .sh`
./core-test.sh $name \
  answers/${name}.txt \
  tested.txt \
  0 0 0 \
  wigs/main.wig wigs/second.wig \
  ../../bwtool paste -skip-NA -consts=3.4,-2.3 main.bw second.bw -o=tested.txt
exit $?
