#!/bin/bash

name=find_max_1
./core-test.sh $name \
  answers/${name}.ans2.bed \
  tested.bed \
  0 0 0 \
  wigs/main.wig \
  ../../bwtool find maxima ../beds/${name}.bed main.bw tested.bed -ave
exit $?
