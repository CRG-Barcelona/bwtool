#ifndef STUFF_H
#define STUFF_H

#ifndef BASICBED_H
#include "basicBed.h"
#endif

/* Misc stuff */

int doubleCountNA(int size, double *array);
/* count the number of NAs in an array */

int doubleWithNASort(int count, double *array);
/* Sort an array of doubles. Return the number of non-NA  */

double doubleWithNAMedianAlreadySorted(int num_non_na, double *array);
/* seems stupid to have but sometimes it's convenient to have the sort and */
/* the median separate */

double doubleWithNAQuantAlreadySorted(int num_non_na, double *array, double quant);
/* this doesn't call Median, which would be the logical thing... but Median */
/* uses bitwise shifting to find the middle array index rather than numerically */
/* **note that this also returns the value at the floored index after mutliplying the */
/*   size by the quantile. For the 0.25 quantile of a size-eight array this will be */
/*   the second value.  Note this isn't quite correct.  Like median, it should */
/*   average the second and third one. */

double doubleWithNAInvQuantAlreadySorted(int num_non_na, double *array, unsigned inv_quant, boolean first);
/* these are special cases that should be more accurate than doubleWithNAQuantAlreadySorted */

double doubleWithNAMedian(int count, double *array);
/* Return median value in array.  This will sort
 * the array as a side effect. */

double doubleWithNAMean(int count, double *array);
/* mean an array assuming they're may be some NA's in there and skip them */

void bedLoadAllReturnFieldCountAndRgbAtLeast3(char *fileName, struct bed **retList, int *retFieldCount, 
					      boolean *retRgb);
/* Load bed of unknown size and return number of fields as well as list of bed items.
 * Ensures that all lines in bed file have same field count.  Also returns whether 
 * column 9 is being used as RGB or not. */
/* This is a copy of bedLoadAllReturnFieldCountAndRgb except it allows bed3 */

#endif /* STUFF_H */
