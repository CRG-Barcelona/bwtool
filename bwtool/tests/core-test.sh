#!/bin/sh

# The core testing script accepts a command like
# 
#    core-test.sh expected-results.txt tested-results.txt decimals wigtype condense A.bw B.bw bwtool paste A.bw B.bw > tested-results.txt
# 
# it will generate bigWigs from wigs, run bwtool, then 
# do a comparison between the two .txt files and return 
# the proper exit code.  

# This is a lousy, flimsy, stupid script but it should work in principle.

name=$1
expected=$2
tested=$3
decimals=$4
wigtype=$5
condense=$6
shift
shift
shift
shift
shift
shift

# expected results file isn't present.  skip this test
if [ ! -e $expected ]; then 
    echo $expected
    exit 77
fi

# we need to avoid race conditions in a parallel check
# e.g. make check -j 4 so make the make a temp dir
tmpdir=`mktemp -d ${name}.XXXX`

bws=""
while [ $1 != "../../bwtool" ]; do
    wig=$1
    bw=`basename $wig .wig`.bw
    sizes=${wig%.wig}.sizes
    # wig isn't present.  skip this test
    if [ ! -e $wig ]; then 
        rm -rf $tmpdir
	exit 77 
    fi
    ./bwmake $sizes $wig ${tmpdir}/$bw
    bws=${bws}" "$bw
    shift
done

cd $tmpdir
$@
# bwtool error'd
if [ $? -gt 0 ]; then
    cd ../
    rm -rf $tmpdir
    exit 2
fi

# output file not generated
if [ ! -e $tested ]; then 
    cd ../
#    rm -rf $tmpdir
    exit 3
fi

# create .wig if needed
if [ $expected != ${expected%.wig} ]; then
   ./wigmake $tested ${tested%.bw}.wig $decimals $wigtype $condense 
   rm -f $tested
   tested=${tested%.bw}.wig
fi 

# expected and tested files not equal 
difference=`diff ../$expected $tested | wc -l`
if [ "$difference" -gt 0 ]; then
    cd ../
    echo "results don't match correct answer"
    mkdir -p fails
    cp ${tmpdir}/$tested fails/${name}-$tested
    rm -fr $tmpdir
    exit 1
fi

cd ../
rm -fr $tmpdir
exit 0
