#!/bin/bash

name=`basename $0 .sh`
./core-test.sh $name \
  answers/${name}.txt \
  tested.txt \
  1 var no \
  wigs/main.wig wigs/second.wig \
  ../../bwtool matrix 5:5 ../beds/agg2.bed main.bw,second.bw tested.txt -starts -keep-bed -decimals=1
exit $?
