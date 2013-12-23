#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "common.h"
#include "obscure.h"
#include "hash.h"
#include "localmem.h"
#include "sqlNum.h"
#include "basicBed.h"
#include "genomeRangeTree.h"
#include "metaBig.h"
#include "bigs.h"
#include "extrema.h"

struct extrema *extrema_new(char *chrom, int chromStart, double val, boolean isMax)
/* create an extrema struct and return it */
{
    struct extrema *ret;
    AllocVar(ret);
    ret->chrom = cloneString(chrom);
    ret->chromStart = chromStart;
    ret->val = val;
    if (isMax)
	ret->min_or_max = '+';
    else
	ret->min_or_max = '-';
    return ret;
}

void extrema_free(struct extrema **pEx)
/* just free one */
{
    if (*pEx)
    {
	struct extrema *ex = *pEx;
	freeMem(ex->chrom);
    }
    freez(pEx);
}

void extrema_free_list(struct extrema **pList)
/* free list */
{
    struct extrema *to_remove;
    while ((to_remove = slPopHead(pList)) != NULL)
    {
	freeMem(to_remove->chrom);
	freez(&to_remove);
    }
    freez(pList);
}

int extrema_cmp(const void *va, const void *vb)
/* Sort first by min_or_max, then by extremity of the value (depends on min_or_max) */
/* i.e. if it's a local max, sort descending, if it's a local min, sort ascending */ 
{
    const struct extrema *a = *((struct extrema **)va);
    const struct extrema *b = *((struct extrema **)vb);
    int dif = 0;
    if (a->min_or_max != b->min_or_max)
	dif = ((a->min_or_max == '+') ? -1 : 1);
    else
    {
	if (((a->min_or_max == '+') && (b->val < a->val)) || ((a->min_or_max == '-') && (b->val > a->val)))
	    dif = -1;
	else if (((a->min_or_max == '+') && (b->val > a->val)) || ((a->min_or_max == '-') && (b->val < a->val)))
	    dif = 1;
    }
    return dif;
}

int extrema_bed_cmp(const void *va, const void *vb)
/* Sort like a bed */
{
    const struct extrema *a = *((struct extrema **)va);
    const struct extrema *b = *((struct extrema **)vb);
    int dif = 0;
    dif = strcmp(a->chrom, b->chrom);
    if (dif != 0)
	return dif;
    return a->chromStart - b->chromStart;
}

static void refine_extrema(struct extrema **pList, unsigned min_sep)
/* sorts by max and value descending, then by min and value ascending */
/* then consecutively adds each extrema depending on whether when extending */
/* the region and adding to a rangeTree, whether that range is already */
/* overlapping something in the tree.  If it does, the extremum is discarded. */
{
    struct extrema *new_list = NULL;
    struct extrema *top;
    int half_sep = min_sep/2;
    slSort(pList, extrema_cmp);
    while (*pList != NULL)
    {
	struct rbTree *rt = rangeTreeNew();
	int start, end;
	char min_or_max;
	while ((top = slPopHead(pList)) != NULL)
	{
	    struct extrema *next = *pList;
	    min_or_max = top->min_or_max;
	    start = top->chromStart - half_sep;
	    end = top->chromStart + half_sep;
	    if (!rangeTreeOverlaps(rt, start, end))
	    {
		rangeTreeAdd(rt, start, end);
		slAddHead(&new_list, top);
	    }
	    else
		extrema_free(&top);
	    if (next && (min_or_max != next->min_or_max))
		break;
	}
	rbTreeFree(&rt);
    }
    *pList = new_list;
    slSort(pList, bedCmp);
}

static void remove_min_or_max(struct extrema **pList, char min_or_max)
/* remove all the extrema in a certain category */
{
    struct extrema *ex;
    struct extrema *new_list = NULL;
    while ((ex = slPopHead(pList)) != NULL)
    {
	if (ex->min_or_max == min_or_max)
	    extrema_free(&ex);
	else
	    slAddHead(&new_list, ex);
    }
    slReverse(&new_list);
    *pList = new_list;
}

static void local_ex(struct extrema **pList, struct perBaseWig *pbw, int i, int j, boolean isMax)
/* add local extrema to list */
{
    struct extrema *newone;
    int middle = i + (j-i)/2;
    newone = extrema_new(pbw->chrom, pbw->chromStart + middle, pbw->data[middle], isMax);
    slAddHead(pList, newone);
}

struct extrema *extrema_find(struct metaBig *mb, int min_sep, enum ex_removal rem)
/* the heart of the algorithm.  does the three-point per-base extrema search */
{
    struct bed *section;
    struct extrema *big_list = NULL;
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbwList = perBaseWigLoadContinue(mb, section->chrom, section->chromStart, 
							      section->chromEnd);
	struct perBaseWig *pbw;
	struct extrema *section_list = NULL;
	for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
	{
	    int size = pbw->chromEnd - pbw->chromStart;
	    if (size > 2)
	    {
		int last = 0;
		int middle = 1;
		int next = 2;
		for (last = 0, middle = 1, next = 2; next < size; next++, middle = next-1, last = middle-1)
		{
		    while ((next+1 < size) && (pbw->data[next] == pbw->data[middle]))
		    {
			middle++;
			next++;
		    }
		    if ((pbw->data[middle] > pbw->data[next]) && (pbw->data[middle] > pbw->data[last]))
			local_ex(&section_list, pbw, last, next, TRUE);
		    else if ((pbw->data[middle] < pbw->data[next]) && (pbw->data[middle] < pbw->data[last]))
			local_ex(&section_list, pbw, last, next, FALSE);
		}
	    }
	    if (min_sep == 0)
		slReverse(&section_list);
	}
	if (min_sep > 0)
	    refine_extrema(&section_list, min_sep);
	big_list = slCat(big_list, section_list);
	perBaseWigFree(&pbwList);
    }
    if (rem == remove_min)
	remove_min_or_max(&big_list, '-');
    else if (rem == remove_max)
	remove_min_or_max(&big_list, '+');
    return big_list;
}

static struct genomeRangeTree *grtFromExtrema(struct extrema *other_list)
/* put the starts in the genomeRangeTree */
{
    struct genomeRangeTree *grt = genomeRangeTreeNew();
    struct extrema *cur;
    for (cur = other_list; cur != NULL; cur = cur->next)
	genomeRangeTreeAdd(grt, cur->chrom, cur->chromStart, cur->chromStart + 1);
    return grt;
}

void extrema_find_shifts(struct extrema *main_list, struct extrema *other_list, unsigned shift)
/* using the two lists and that shift value, perform nearest-neighbor search of extrema */
/* in other list.  the result is to just change the main_list's values from value to distance. */
/* -1 distance indicates other extremum not seen in range */
{
    struct genomeRangeTree *grt = grtFromExtrema(other_list);
    struct extrema *cur;
    for (cur = main_list; cur != NULL; cur = cur->next)
    {
	int s = cur->chromStart - shift;
	int e = cur->chromStart + shift;
	if (s < 0)
	    s = 0;
	struct range *list = genomeRangeTreeAllOverlapping(grt, cur->chrom, s, e);
	struct range *rn;
	int lowest = shift + 1;
	for (rn = list; rn != NULL; rn = rn->next)
	{
	    int dist = abs(rn->start - cur->chromStart);
	    if (dist < lowest)
		lowest = dist;
	}
	if (lowest < shift)
	    cur->val = lowest;
	else
	    cur->val = -1;
    }
}
