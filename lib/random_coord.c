/* module to deal with random coordinates across chromosomes */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "common.h"
#include "obscure.h"
#include "random_coord.h"
#include "bigs.h"

#define NANUM sqrt(-1)

#ifdef HAVE_LIBGSL
#include <gsl/gsl_rng.h>
#endif

/* "private" struct declaration */
struct cat_chrom
{
    struct cat_chrom *next;
    char *chrom;
    unsigned chromStart;
    unsigned chromEnd;
    unsigned size;
    unsigned long start;
    unsigned long end;
};

static boolean convert_bed_to_range(struct random_coord *rc, struct bed *bed, int *start, int *end)
/* convert the given bed to the catted-space range coordinates, returning the start */
/* and end inside pointers in parameter list. if chromosome is not one of the ones in */
/* the random_coord, return FALSE so it may be skipped */
{    
    struct cat_chrom *cc = hashFindVal(rc->lookup_by_chrom, bed->chrom);
    if (!cc)
	return FALSE;
    *start = cc->start + bed->chromStart;
    *end = *start + (bed->chromEnd - bed->chromStart);
    return TRUE;
}

void random_coord_add_black(struct random_coord *rc, struct bed *blacklist)
/* merge more bad regions with what's already in the current blacklist */
{
    struct bed *bed;
    for (bed = blacklist; bed != NULL; bed = bed->next)
    {
	int start, end;
	if (convert_bed_to_range(rc, bed, &start, &end))
	    rangeTreeAdd(rc->blacklist, start, end);
    }
}

struct random_coord *random_coord_init(struct hash *chrom_sizes, struct bed *blacklist)
/* Initialize struct using a typical chrom_sizes hash and an optional bedfile that's */
/* already been read in. The hash is cloned, so the original needs free'ing separately */
{
    struct random_coord *ret;
    struct hashEl *list = hashElListHash(chrom_sizes), *el;
    struct cat_chrom *cc;
    AllocVar(ret);
    ret->lookup_by_coord = rangeTreeNew();
    ret->lookup_by_chrom = hashNew(chrom_sizes->powerOfTwoSize);
    ret->length = 0;
    ret->blacklist = rangeTreeNew();
    /* build initial list and count up length */
    for (el = list; el != NULL; el = el->next)
    {
	struct cat_chrom *new_cc;
	AllocVar(new_cc);
	new_cc->chrom = cloneString(el->name);
	new_cc->size = ptToInt(el->val);
	new_cc->start = ret->length;	
	ret->length += new_cc->size;
	new_cc->end = ret->length;
	hashAdd(ret->lookup_by_chrom, new_cc->chrom, new_cc);
	slAddHead(&ret->chroms, new_cc);
    }
    slReverse(&ret->chroms);
    /* make reverse lookup list */
    for (cc = ret->chroms; cc != NULL; cc = cc->next)
	rangeTreeAddVal(ret->lookup_by_coord, cc->start, cc->end, cc, NULL);
    random_coord_add_black(ret, blacklist);
    hashElFreeList(&list);
    return ret;
}

static void cat_chrom_free_list(struct cat_chrom **pList)
/* free up the list */
{
    struct cat_chrom *el;
    while ((el = slPopHead(pList)) != NULL)
    {
	freeMem(el->chrom);
	freez(&el);
    }
}

void random_coord_free(struct random_coord **pRc)
/* free-up the random coordinate struct */
{
    struct random_coord *rc = *pRc;
    if (rc)
    {
	cat_chrom_free_list(&rc->chroms);
	rangeTreeFree(&rc->lookup_by_coord);
	rangeTreeFree(&rc->blacklist);
	hashFree(&rc->lookup_by_chrom);
	freez(pRc);
    }
}

struct bed *random_bed(struct random_coord *rc, int size, int rand_coord)
/* retrieve a region of a given size. random number should be generated before-hand.  */
/* returns NULL if random coordinate is for some reason out of range, or if the */
/* rsulting bed is out of range, encloses a boundary, or overlaps a blacklisted range */
{
    struct cat_chrom *cc;
    struct bed *bed;
    struct range *lookedup = NULL;
    if (((unsigned long)rand_coord + size > rc->length) || (rangeTreeOverlaps(rc->blacklist, rand_coord, rand_coord + size)))
	return NULL;
    lookedup = rangeTreeFindEnclosing(rc->lookup_by_coord, rand_coord, rand_coord + size);
    if (!lookedup)
	return NULL;
    cc = (struct cat_chrom *)lookedup->val;
    AllocVar(bed);
    bed->chrom = cloneString(cc->chrom);
    bed->chromStart = rand_coord - cc->start;
    bed->chromEnd = bed->chromStart + size;
    return bed;
}

#ifdef HAVE_LIBGSL

struct perBaseWig *random_pbw_list(int size, int N, struct metaBig *mb, double NA_perc, 
				   double fill, struct bed *blacklist, unsigned seed)
/* retrieve a list random regions' data from a bigWig */
{
    struct perBaseWig *pbwList = NULL;
    const gsl_rng_type *T;
    gsl_rng *r;
    struct random_coord *rc = NULL;
    unsigned long max = 0;
    unsigned kept = 0;
    double na = NANUM;
    unsigned long int s = time(0);
    rc = random_coord_init(mb->chromSizeHash, blacklist);
    max = rc->length - size;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    if (seed != 0)
	gsl_rng_set(r,(unsigned long)seed);
    else 
	gsl_rng_set(r,s);
    while (kept < N)
    {
	struct bed *bed = NULL;
	while (bed == NULL)
	{
	    unsigned long rand = gsl_rng_uniform_int(r, max);
	    bed = random_bed(rc, size, rand);
	}
	struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, bed->chrom, bed->chromStart, bed->chromEnd, FALSE, fill);
	if (isnan(fill))
	{
	    int size_pbw = pbw->len;
	    int count_NA = 0;
	    int i;
	    for (i = 0; i < size_pbw; i++)
	    {
		if (isnan(pbw->data[i]))
		    count_NA++;
	    }
	    if ((double)count_NA/size_pbw <= NA_perc)
	    {
		slAddHead(&pbwList, pbw);
		kept++;
	    }
	}
	else 
	{
	    slAddHead(&pbwList, pbw);
	    kept++;
	}
	bedFree(&bed);
    }
    gsl_rng_free(r);
    random_coord_free(&rc);
    return pbwList;
}

struct bed *random_bed_list(int size, int N, struct metaBig *mb, double NA_perc, 
			    double NA, struct bed *blacklist)
/* same thing as the pbw list, but return only coordinates */
{
    double na = NANUM;
    struct perBaseWig *pbwList = random_pbw_list(size, N, mb, NA_perc, na, blacklist, 0);
    struct perBaseWig *pbw;
    struct bed *bedList = NULL;
    for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
    {
	struct bed *bed;
	AllocVar(bed);
	bed->chrom = cloneString(pbw->chrom);
	bed->chromStart = pbw->chromStart;
	bed->chromEnd = pbw->chromEnd;
	slAddHead(&bedList, bed);
    }
    slReverse(&bedList);
    perBaseWigFreeList(&pbwList);
    return bedList;
}

#else
/* not HAVE_LIBGSL */

struct perBaseWig *random_pbw_list(int size, int N, struct metaBig *mb, double NA_perc, 
				   double fill, struct bed *blacklist, unsigned seed)
{
    return NULL;
}

struct bed *random_bed_list(int size, int N, struct metaBig *mb, double NA_perc, 
			    double NA, struct bed *blacklist)
{
    return NULL;
}

#endif
/* HAVE_LIBGSL */
