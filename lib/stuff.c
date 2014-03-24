#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "common.h"
#include "stuff.h"

#define NANUM sqrt(-1)

static int doubleWithNACmp(const void *va, const void *vb)
/* Compare function to sort array of doubles. */
{
    const double *a = va;
    const double *b = vb;
    double diff;
    if (isnan(*a) && isnan(*b))
	return 0;
    if (isnan(*a))
	return 1;
    if (isnan(*b))
	return -1;
    diff = *a - *b;
    if (diff < 0)
	return -1;
    if (diff > 0)
	return 1;
    return 0;
}

int doubleWithNASort(int count, double *array)
/* Sort an array of doubles. Return the number of non-NA  */
{
    int size = 0;
    if (count > 1)
	qsort(array, count, sizeof(array[0]), doubleWithNACmp);
    else
	return count;
    while ((size < count) && (!isnan(array[size])))
	size++;
    return size;
}

double doubleWithNAMedianAlreadySorted(int num_non_na, double *array)
/* seems stupid to have but sometimes it's convenient to have the sort and */
/* the median separate */
{
    double median;
    if (num_non_na == 0)
	median = NANUM;
    else if ((num_non_na&1) == 1)
	median = array[num_non_na>>1];
    else
    {
	num_non_na >>= 1;
	median = (array[num_non_na] + array[num_non_na-1]) * 0.5;
    }
    return median;
}

double doubleWithNAQuantAlreadySorted(int num_non_na, double *array, double quant)
/* this doesn't call Median, which would be the logical thing... but Median */
/* uses bitwise shifting to find the middle array index rather than numerically */
/* **note that this also returns the value at the floored index after mutliplying the */
/*   size by the quantile. For the 0.25 quantile of a size-eight array this will be */
/*   the second value.  Note this isn't quite correct.  Like median, it should */
/*   average the second and third one. */
{
    int ix = quant * num_non_na;
    if ((quant <= 0) || (quant >= 1))
	errAbort("bad quantile specified:  pick between 0-1.0");
    if (ix == num_non_na)
	ix--;
    return array[ix];
}

double doubleWithNAInvQuantAlreadySorted(int num_non_na, double *array, unsigned inv_quant, boolean first)
/* these are special cases that should be more accurate than doubleWithNAQuantAlreadySorted */
{
    if (inv_quant < 2)
	errAbort("Need to specify larger than 2 for inverse quantile");
    float quant = (first) ? (1/(float)inv_quant) : 1-1/(float)inv_quant;
    int ix = (int)(num_non_na * quant);
    if (num_non_na % inv_quant == 0)
	return (array[ix-1] + array[ix])/2;
    return array[ix];
}

double doubleWithNAMedian(int count, double *array)
/* Return median value in array.  This will sort
 * the array as a side effect. */
{
    int size = doubleWithNASort(count, array);
    return doubleWithNAMedianAlreadySorted(size, array);
}

double doubleWithNAMean(int count, double *array)
/* mean an array assuming they're may be some NA's in there and skip them */
{
    long double sum = 0;
    int i;
    int num_non_na = 0;
    for (i = 0; i < count; i++)
    {
	if (!isnan(array[i]))
	{
	    num_non_na++;
	    sum += array[i];
	}
    }
    if (num_non_na == 0)
	return NANUM;
    return (double)(sum/num_non_na);
}

void bedLoadAllReturnFieldCountAndRgbAtLeast3(char *fileName, struct bed **retList, int *retFieldCount, 
    boolean *retRgb)
/* Load bed of unknown size and return number of fields as well as list of bed items.
 * Ensures that all lines in bed file have same field count.  Also returns whether 
 * column 9 is being used as RGB or not. */
/* This is a copy of bedLoadAllReturnFieldCountAndRgb except it allows bed3 */
{
struct bed *list = NULL;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *line, *row[bedKnownFields];
int fieldCount = 0;
boolean isRgb = FALSE;

while (lineFileNextReal(lf, &line))
    {
    int numFields = chopByWhite(line, row, ArraySize(row));
    if (numFields < 3)
	errAbort("file %s doesn't appear to be in bed format. At least 3 fields required, got %d", 
		fileName, numFields);
    if (fieldCount == 0)
	{
        fieldCount = numFields;
	if (fieldCount >= 9)
	    isRgb =  (strchr(row[8], ',') != NULL);
	}
    else
        if (fieldCount != numFields)
	    errAbort("Inconsistent number of fields in file. %d on line %d of %s, %d previously.",
	        numFields, lf->lineIx, lf->fileName, fieldCount);
    slAddHead(&list, bedLoadN(row, fieldCount));
    }
lineFileClose(&lf);
slReverse(&list);
*retList = list;
*retFieldCount = fieldCount;
if (retRgb != NULL)
   *retRgb = isRgb;
}
