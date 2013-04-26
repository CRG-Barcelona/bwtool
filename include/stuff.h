#ifndef STUFF_H
#define STUFF_H

/* Misc stuff */

int doubleCountNA(int size, double *array);
/* count the number of NAs in an array */

int doubleWithNASort(int count, double *array);
/* Sort an array of doubles. Return the number of non-NA  */

double doubleWithNAMedianAlreadySorted(int num_non_na, double *array);
/* seems stupid to have but sometimes it's convenient to have the sort and */
/* the median separate */

double doubleWithNAMedian(int count, double *array);
/* Return median value in array.  This will sort
 * the array as a side effect. */

double doubleWithNAMean(int count, double *array);
/* mean an array assuming they're may be some NA's in there and skip them */

#endif /* STUFF_H */
