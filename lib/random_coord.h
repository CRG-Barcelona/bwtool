/* module to deal with random coordinates across chromosomes */

#ifndef RANDOM_COORD_H
#define RANDOM_COORD_H

#ifndef HASH_H
#include "hash.h"
#endif

#ifndef BASICBED_H
#include "basicBed.h"
#endif

#ifndef RANGETREE_H
#include "rangeTree.h"
#endif

#ifndef METABIG_H
#include "metaBig.h"
#endif

typedef struct { } cat_chrom;
/* declared "privately" within the .c */

struct random_coord
/* hold all the info needed */
{
    struct cat_chrom *chroms;
    unsigned long length;
    struct rbTree *blacklist;        /* regions in catted-space not to be overlapped */
    struct hash *lookup_by_chrom;    /* given a chromosome name, find the associated start */
    struct rbTree *lookup_by_coord;  /* used to find chroms from catted-space */
};

int random_in_range(unsigned long min, unsigned long max);
/* Would like a semi-open interval [min, max) */

struct random_coord *random_coord_init(struct hash *chrom_sizes, struct bed *blacklist);
/* Initialize struct using a typical chrom_sizes hash and an optional bedfile that's */
/* already been read in. The hash is cloned, so the original needs free'ing separately */

void random_coord_free(struct random_coord **pRc);
/* free-up the random coordinate struct */

void random_coord_add_black(struct random_coord *rc, struct bed *blacklist);
/* merge more bad regions with what's already in the current blacklist */

struct bed *random_bed(struct random_coord *rc, int size, int rand_coord);
/* retrieve a region of a given size. random number should be generated before-hand.  */
/* returns NULL if random coordinate is for some reason out of range, or if the */
/* rsulting bed is out of range, encloses a boundary, or overlaps a blacklisted range */

struct perBaseWig *random_pbw_list(int size, int N, struct metaBig *mb, double NA_perc, 
				   double fill, struct bed *blacklist, unsigned seed);
/* retrieve a list random regions' data from a bigBed */

struct bed *random_bed_list(int size, int N, struct metaBig *mb, double NA_perc, 
			    double NA, struct bed *blacklist);
/* same thing as the pbw list, but return only coordinates */

#endif /* RANDOM_COORD_H */
