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
