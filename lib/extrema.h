#ifndef EXTREMA_H
#define EXTREMA_H

enum ex_removal
/* Just an enum to track the -minima and -maxima options */
{
    no_removal = 0,
    remove_min = 1,
    remove_max = 2,
};

struct extrema 
/* same thing as a bed6, just internally name is a double in order */
/* to skip nonsense string conversions later */
{
    struct extrema *next;
    char *chrom;
    int chromStart;
    double val;
    char min_or_max;
};

struct extrema *new_ex(char *chrom, int chromStart, double val, boolean isMax);
/* create an extrema struct and return it */

void extrema_free(struct extrema **pEx);
/* just free one */

void extrema_free_list(struct extrema **pList);
/* free list */

int extrema_cmp(const void *va, const void *vb);
/* Sort first by min_or_max, then by extremity of the value (depends on min_or_max) */
/* i.e. if it's a local max, sort descending, if it's a local min, sort ascending */ 

int extrema_bed_cmp(const void *va, const void *vb);
/* Sort first by min_or_max, then by extremity of the value (depends on min_or_max) */
/* i.e. if it's a local max, sort descending, if it's a local min, sort ascending */ 

struct extrema *extrema_find(struct metaBig *mb, int min_sep, enum ex_removal rem);
/* the heart of the algorithm.  does the three-point per-base extrema search */

void extrema_find_shifts(struct extrema *main_list, struct extrema *other_list, unsigned shift);
/* using the two lists and that shift value, perform nearest-neighbor search of extrema */
/* in other list.  the result is to just change the main_list's values from value to distance. */
/* -1 distance indicates other extremum not seen in range */

#endif /* define EXTREMA_H */
