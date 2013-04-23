/* bwtool_distance - distance between signals  */

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "basicBed.h"
#include "bigWig.h"
#include "bigs.h"
#include "bwtool.h"

void usage_distance()
/* Explain usage of distance calculation program and exit. */
{
errAbort(
    "bwtool distance - calculate distance between signals\n")
    "usage:\n"
    "   bwtool distance type input.bw[:chr:start-end] output.txt\n" 
    "where \"type\" is \"pearson-corr\""
    "options:\n"
    "   -window=n      calculate distance over sliding windows\n"
  );
}


