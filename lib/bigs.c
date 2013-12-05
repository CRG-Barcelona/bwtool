#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "common.h"
#include "obscure.h"
#include "hash.h"
#include "linefile.h"
#include "localmem.h"
#include "sqlNum.h"
#include "sig.h"
#include "dystring.h"
#include "basicBed.h"
#include "bigBed.h"
#include "bigWig.h"
#include "hmmstats.h"
#include "bigs.h"

#include <limits.h>
#include <math.h>

#define NANUM sqrt(-1)

enum wigOutType get_wig_out_type(char *option)
/* scan the option and return it */
{
    char *wig_type = NULL;
    if (!option)
	wig_type = "fix";
    else 
	wig_type = option;
    if (sameString("bg", wig_type))
	return bedGraphOut;
    else if (sameString("fix", wig_type))
	return fixStepOut;
    else if (sameString("var", wig_type))
	return varStepOut; 
    return invalidOut;
}

double bigWigMean(struct bbiFile *bw)
/* return the mean value of a bigWig */
{
    double na = NANUM;
    struct bbiSummaryElement bbs = bbiTotalSummary(bw);
    if (bbs.validCount == 0)
	return na;
    else
	return bbs.sumData/bbs.validCount;
}

double bigWigStd(struct bbiFile *bw)
/* return the mean value of a bigWig */
{
    double na = NANUM;
    struct bbiSummaryElement bbs = bbiTotalSummary(bw);
    if (bbs.validCount == 0)
	return na;
    else
	return calcStdFromSums(bbs.sumData, bbs.sumSquares, bbs.validCount);
}

int perBaseWigLabelCmp(const void *a, const void *b)
/* for sorting after clustering */
{
    const struct perBaseWig *pbw_a = *((struct perBaseWig **)a);
    const struct perBaseWig *pbw_b = *((struct perBaseWig **)b);    
    int diff = pbw_a->label - pbw_b->label;
    if (diff == 0)
    {
	double dist_diff = pbw_b->cent_distance - pbw_a->cent_distance;
	if (dist_diff < 0)
	    diff = 1;
	else if (dist_diff > 0)
	    diff = -1;
    }
    return diff;
}

struct perBaseWig *alloc_perBaseWig(char *chrom, int start, int end)
/* simply allocate the perBaseWig. this is filled with NA values */
/* it's best to call this one if we are going to read a bigWig */
{
    struct perBaseWig *pbw;
    const double na = NANUM;
    int size = end - start;
    int i;
    AllocVar(pbw);
    pbw->chrom = cloneString(chrom);
    pbw->chromStart = start;
    pbw->chromEnd = end;
    pbw->len = size;
    if (size > 0)
    {
	AllocArray(pbw->data, size);
	for (i = 0; i < size; i++)
	    pbw->data[i] = na;
    }
    return pbw;
}

struct perBaseWig *alloc_zero_perBaseWig(char *chrom, int start, int end)
/* simply allocate the perBaseWig. this is filled with zeros */
/* it may be best to call this one before writing a wig */
{
    struct perBaseWig *pbw = alloc_perBaseWig(chrom, start, end);
    int i;
    for (i = 0; i < pbw->len; i++)
	pbw->data[i] = 0;
    return pbw;
}

struct perBaseWig *alloc_fill_perBaseWig(char *chrom, int start, int end, double fill)
/* fill the pbw with a given value instead of NA */
{
    struct perBaseWig *pbw = alloc_perBaseWig(chrom, start, end);
    if (!isnan(fill))
    {
	int i;
	for (i = 0; i < pbw->len; i++)
	    pbw->data[i] = fill;
    }
    return pbw;    
}

struct perBaseWig *perBaseWigClone(struct perBaseWig *pbw)
/* Clone the structure exactly */
{
    struct perBaseWig *clone;
    int size = pbw->chromEnd - pbw->chromStart;
    int i;
    struct bed *pbw_bed;
    AllocVar(clone);
    clone->chrom = cloneString(pbw->chrom);
    clone->chromStart = pbw->chromStart;
    clone->chromEnd = pbw->chromEnd;
    clone->len = size;
    AllocArray(clone->data, size);
    for (i = 0; i < size; i++)
	clone->data[i] = pbw->data[i];
    clone->subsections = NULL;
    for (pbw_bed = pbw->subsections; pbw_bed != NULL; pbw_bed = pbw_bed->next)
    {
	struct bed *newbed;
	AllocVar(newbed);
	newbed->chrom = cloneString(pbw_bed->chrom);
	newbed->chromStart = pbw_bed->chromStart;
	newbed->chromEnd = pbw_bed->chromEnd;
	slAddHead(&clone->subsections, newbed);
    }
    slReverse(&clone->subsections);
    return clone;
}

static char *seq_name_disassemble(struct dnaSeq *seq, unsigned *chromStart, unsigned *chromEnd)
    /* This chops a sequence name up at the first ':' and fills the chromStart/chromEnd */
    /* with the expected values */
{
    unsigned start = 0;
    unsigned end = seq->size;
    char *s = strchr(seq->name, ':');
    if (s)
    {
	*s++ = '\0';
	char *e = strchr(s, '-');
	if (e)
	{
	    *e++ = '\0';
	    end = sqlUnsigned(e);
	}
	start = sqlUnsigned(s);
    }
    *chromStart = (int)start;
    *chromEnd = (int)end;
    return cloneString(seq->name);
}

static struct bed *seq_subsection_list(struct dnaSeq *seq, boolean skipN)
    /* find all the good sections of sequence */
{
    struct bed *list = NULL;
    if (skipN)
    {
	int start = 0;
	int end = 0;
	while (end < seq->size)
	{
	    while ((end < seq->size) && (seq->dna[end] == 'N'))
		end++;
	    start = end;
	    while ((end < seq->size) && (seq->dna[end] != 'N'))
		end++;
	    if (end - start > 0)
	    {
		struct bed *section;
		AllocVar(section);
		section->chrom = cloneString(seq->name);
		section->chromStart = start;
		section->chromEnd = end;
		slAddHead(&list, section);
	    }
	}
	slReverse(&list);
    }
    else
    {
	struct bed *wholeseq;
	AllocVar(wholeseq);
	wholeseq->chrom = cloneString(seq->name);
	wholeseq->chromStart = 0;
	wholeseq->chromEnd = seq->size;
	slAddHead(&list, wholeseq);
    }
    return list;
}

struct perBaseWig *alloc_perBaseWig_matchingSequence(struct dnaSeq *seq, boolean skipN)
/* allocate a perBaseWig to match the length of the dnaSeq.  Optionally choose to skip */
/* N bases by making a subsections bed list that avoids them. One feature this */
/* function has is that the name of the sequence is something like "chrom:start-end" */
/* in 0-based coordinates. If so, the chromStart/chromEnd are set to match the coordinates in */
/* the name.  If not, the chromStart will be 0, and the chromEnd will be seq->size. */
/* (this is used by symcurv) */
{
    struct perBaseWig *pbw;
    int size = seq->size;
    unsigned chromStart;
    unsigned chromEnd;
    char *chrom = seq_name_disassemble(seq, &chromStart, &chromEnd);
    pbw = alloc_perBaseWig(chrom, (int)chromStart, (int)chromEnd);
    pbw->subsections = seq_subsection_list(seq, skipN);
    AllocArray(pbw->data, size);
    return pbw;
}

void perBaseWigFree(struct perBaseWig **pRegion)
/* Free-up a perBaseWig */
{
    struct perBaseWig *pbw = *pRegion;
    if (!pRegion || !pbw)
	return;
    if (pbw->subsections)
	bedFreeList(&pbw->subsections);
    if (pbw->name)
	freeMem(pbw->name);
    freeMem(pbw->chrom);
    freez(&pbw->data);
    freez(pRegion);
}

void perBaseWigFreeList(struct perBaseWig **pList)
/* Free up a list of perBaseWigs */
{
    struct perBaseWig *reg;
    while ((reg = slPopHead(pList)) != NULL)
	perBaseWigFree(&reg);
    freez(pList);
}

struct perBaseWig *perBaseWigLoadContinue(struct metaBig *mb, char *chrom, int start, int end)
/* load a perBaseWig from a wig/bigWig that's already open */
{
    if (mb->type != isaBigWig)
	return NULL;
    struct perBaseWig *list = NULL;
    struct lm *lm = lmInit(0);
    struct bbiInterval *intervals = bigWigIntervalQuery(mb->big.bbi, chrom, start, end, lm);
    struct bbiInterval *bbStart = intervals, *bbEnd;
    while (bbStart != NULL)
    {
	struct perBaseWig *region;
	struct bbiInterval *cur;
	int size = 0;
	int i = 0;
	bbEnd = bbStart;
	/* loop until discontinuity detected */
	while ((bbEnd->next != NULL) && (bbEnd->end == bbEnd->next->start))
	    bbEnd = bbEnd->next;
	region = alloc_perBaseWig(chrom, bbStart->start, bbEnd->end);
	for (cur = bbStart; cur != bbEnd->next; cur = cur->next)
	{
	    int j;
	    for (j = cur->start; j < cur->end; j++)
		region->data[i++] = cur->val;
	}
	slAddHead(&list, region);
	bbStart = bbEnd->next;
    }
    lmCleanup(&lm);
    slReverse(&list);
    return list;
}

struct perBaseWig *perBaseWigLoad(char *wigFile, char *chrom, int start, int end)
/* Load all the regions from a wig or bigWig into a list of arrays basically. */
{
    struct metaBig *mb = metaBigOpen(wigFile, NULL);
    if (mb->type != isaBigWig)
    {
	metaBigClose(&mb);
	return NULL;
    }
    struct perBaseWig *list;
    list = perBaseWigLoadContinue(mb, chrom, start, end);
    metaBigClose(&mb);
    return list;
}

static void chromOob(struct metaBig *mb, char *chrom, int *start, int *end)
{
    int csize = hashIntVal(mb->chromSizeHash, chrom);
    if (*start < 0)
	*start = 0;
    if (*end > csize)
	*end = csize;
}

static void perBaseWigLoadHugeContinue(struct metaBig *mb, struct perBaseWig *big_pbw, int *big_offset, struct bed *section)
/* Load all the regions into one perBaseWig, but with gaps filled  */
/* in with NA value */
{
    struct perBaseWig *list = NULL;
    double na = NANUM;
    int s = section->chromStart;
    int e = section->chromEnd;
    if (!hashFindVal(mb->chromSizeHash, section->chrom))
    {
	/* if the chrom isn't in the bigWig's chrom-size hash, skip over */
	*big_offset += e - s;
    }
    else
    {
	struct perBaseWig *pbw;
	chromOob(mb, section->chrom, &s, &e);
	list = perBaseWigLoadContinue(mb, section->chrom, s, e);
	if (list)
	{
	    for (pbw = list; pbw != NULL; pbw = pbw->next)
	    {
		int j;
		for (j = 0; j < pbw->len; j++)
		    big_pbw->data[*big_offset + j] = pbw->data[j];
		*big_offset += pbw->len;
	    }
	    perBaseWigFreeList(&list);
	}
    }
}

struct perBaseWig *perBaseWigLoadHuge(struct metaBig *mb, struct bed *regions)
/* Load a huge pbw, gaps removed */
{
    long supposed_size = 0;
    struct bed *bed;
    struct perBaseWig *big_pbw = NULL; 
    int big_offset = 0;
    for (bed = regions; bed != NULL; bed = bed->next)
	supposed_size += bed->chromEnd - bed->chromStart;
    if (supposed_size > powl(2,31))
	errAbort("Requested regions sum to greater than 2^31 = 2,147,483,648 bases. The current implementation is restricted to fewer than this");
    big_pbw = alloc_perBaseWig("various", 0, supposed_size);
    for (bed = regions; bed != NULL; bed = bed->next)
	perBaseWigLoadHugeContinue(mb, big_pbw, &big_offset, bed);
    big_pbw->total_coverage = (unsigned)big_offset;
    return big_pbw;
}

struct perBaseWig *perBaseWigLoadSingleContinue(struct metaBig *mb, char *chrom, 
						int start, int end, boolean reverse, double fill)
/* Load all the regions into one perBaseWig, but with gaps filled  */
/* in with NA value */
{
    if (mb->type != isaBigWig)
	errAbort("tried to load data from a non-bigWig file");
    struct perBaseWig *list;
    struct perBaseWig *region;
    struct perBaseWig *wholething = NULL;
    int size = end - start;
    int i, j;
    int s = start, e = end;
    double na = NANUM;
    if (!hashFindVal(mb->chromSizeHash, chrom))
    {
	/* if the chrom isn't in the bigWig's chrom-size hash, return values of NA */
	wholething = alloc_fill_perBaseWig(chrom, start, end, fill);
	return wholething;
    }
    chromOob(mb, chrom, &s, &e);
    list = perBaseWigLoadContinue(mb, chrom, s, e);
    wholething = alloc_fill_perBaseWig(chrom, start, end, fill);
    if (list)
    {
	for (region = list; region != NULL; region = region->next)
	{
	    int offset = region->chromStart - wholething->chromStart;
	    for (j = 0; j < region->chromEnd - region->chromStart; j++)
		wholething->data[offset + j] = region->data[j];
	}
	perBaseWigFreeList(&list);
    }
    if (reverse)
    {
	double swap;
	for (i = 0; i < (size/2); i++)
	{
	    j = (size-1)-i;
	    swap = wholething->data[i];
	    wholething->data[i] = wholething->data[j];
	    wholething->data[j] = swap;
	}
    }
    return wholething;
}

struct perBaseWig *perBaseWigLoadSingleMetaContinue(struct metaBig *mb, char *chrom, 
						int start, int end, boolean reverse, double fill, int size)
/* Load all the regions into one perBaseWig, but with gaps filled  */
/* in with NA value */
{
    struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, chrom, start, end, reverse, fill);
    struct perBaseWig *meta = alloc_fill_perBaseWig(chrom, start, start+size, fill);
    float step = (float)pbw->len/meta->len;
    int i;
    for (i = 0; i < meta->len; i++)
    {
	float pbw_s = step * i;
	float pbw_e = pbw_s + step;
	if ((int)pbw_e - (int)pbw_s == 0)
	    meta->data[i] = pbw->data[(int)pbw_s];
	else
	{
	    double area = 0;
	    double len = 0;
	    int j;
	    int firstix = floor(pbw_s);
	    int between_s = ceil(pbw_s);
	    int between_e = floor(pbw_e);
	    int lastixp1 = ceil(pbw_e);
	    if ((between_s > firstix) && (!isnan(pbw->data[firstix])))
	    {
		area += ((float)between_s - pbw_s) * pbw->data[firstix];
		len += (float)between_s - pbw_s;
	    }
	    for (j = between_s; j < between_e; j++)
		if (!isnan(pbw->data[j]))
		{
		    area += pbw->data[j];
		    len += 1.0;
		}
	    if ((between_e < pbw->len) && (lastixp1 > between_e) && (!isnan(pbw->data[between_e])))
	    {
		area += (pbw_e - (float)between_e) * pbw->data[between_e];
		len += (pbw_e - (float)between_e);
	    }
	    if (len > 0)
		meta->data[i] = area/len;
	    else meta->data[i] = fill;
	}
    }
    perBaseWigFree(&pbw);
    return meta;
}

struct perBaseWig *perBaseWigLoadSingle(char *wigFile, char *chrom, int start, int end, boolean reverse, double fill)
/* Load all the regions into one perBaseWig, but with gaps filled  */
/* in with NA value */
{
    struct metaBig *mb = metaBigOpen(wigFile, NULL);
    if (mb->type != isaBigWig)
    {
	metaBigClose(&mb);
	return NULL;
    }
    struct perBaseWig *list;
    list = perBaseWigLoadSingleContinue(mb, chrom, start, end, reverse, fill);
    metaBigClose(&mb);
    return list;
}

struct perBaseMatrix *load_perBaseMatrix(struct metaBig *mb, struct bed6 *regions, double fill)
/* loading a wig matrix from a metaBig with regions all the same size.  It should be noted */
/* that the regions may include negative and out-of-bounds coordinates.  Out-of-bounds data is */
/* NA in the matrix. */ 
{
    struct perBaseMatrix *pmat;
    struct bed6 *item;
    int i;
    if (!mb->sections)
	return NULL;
    item = regions;
    if (!item)
	return NULL;
    AllocVar(pmat);
    pmat->nrow = slCount(regions);
    pmat->ncol = item->chromEnd - item->chromStart;
    AllocArray(pmat->array, pmat->nrow);
    AllocArray(pmat->matrix, pmat->nrow);
    for (i = 0; (i < pmat->nrow) && (item != NULL); i++, item = item->next)
    {
	struct perBaseWig *pbw;
	pbw = perBaseWigLoadSingleContinue(mb, item->chrom, item->chromStart, item->chromEnd, item->strand[0]=='-', fill);
	pbw->name = cloneString(item->name);
	pbw->score = item->score;
	pbw->len = pmat->ncol;
	pbw->strand[0] = item->strand[0];
	pmat->array[i] = pbw;
	pmat->matrix[i] = pbw->data;
    }
    return pmat;
}

struct perBaseMatrix *load_ave_perBaseMatrix(struct metaBig *mb, struct bed6 *regions, int tile_size, double fill)
/* the matrix is tiled averages instead the values at each base */
{
    struct perBaseMatrix *pmat;
    struct bed6 *item;
    int i, j, k;
    double na = NANUM;
    if (!mb->sections)
	return NULL;
    item = regions;
    if (!regions || (((regions->chromEnd - regions->chromStart) % tile_size) != 0))
	return NULL;
    AllocVar(pmat);
    pmat->nrow = slCount(regions);
    pmat->ncol = item->chromEnd - item->chromStart;
    pmat->ncol /= tile_size;
    AllocArray(pmat->array, pmat->nrow);
    AllocArray(pmat->matrix, pmat->nrow);
    for (i = 0; (i < pmat->nrow) && (item != NULL); i++, item = item->next)
    {
	struct perBaseWig *pbw;
	struct perBaseWig *small = alloc_fill_perBaseWig(item->chrom, item->chromStart, item->chromStart + pmat->ncol, fill);
	small->chromEnd = item->chromEnd;
	pbw = perBaseWigLoadSingleContinue(mb, item->chrom, item->chromStart, item->chromEnd, item->strand[0]=='-', fill);
	for (j = 0; j < pmat->ncol; j++)
	{
	    double sum = 0, ave;
	    int n = 0;
	    for (k = j*tile_size; k < j*tile_size+tile_size; k++)
		if (!isnan(pbw->data[k]))
		{
		    sum += pbw->data[k];
		    n++;
		}
	    if (n > 0)
		ave = sum / n;
	    else
		ave = na;
	    small->data[j] = ave;
	}
	small->name = cloneString(item->name);
	small->score = item->score;
	small->strand[0] = item->strand[0];
	pmat->array[i] = small;
	pmat->matrix[i] = small->data;
	perBaseWigFree(&pbw);
    }
    return pmat;    
}

struct perBaseMatrix *load_meta_perBaseMatrix(struct metaBig *mb, struct bed6 *regions, int meta_size, double fill)
/* Load variably-sized regions into fixed-sized matrix */
{
    struct perBaseMatrix *pmat;
    struct bed6 *item;
    int i;
    if (!mb->sections)
	return NULL;
    item = regions;
    if (!item)
	return NULL;
    AllocVar(pmat);
    pmat->nrow = slCount(regions);
    pmat->ncol = meta_size;
    AllocArray(pmat->array, pmat->nrow);
    AllocArray(pmat->matrix, pmat->nrow);
    for (i = 0; (i < pmat->nrow) && (item != NULL); i++, item = item->next)
    {
	struct perBaseWig *pbw;
	pbw = perBaseWigLoadSingleMetaContinue(mb, item->chrom, item->chromStart, item->chromEnd, item->strand[0]=='-', fill, meta_size);
	pbw->name = cloneString(item->name);
	pbw->score = item->score;
	pbw->len = pmat->ncol;
	pbw->strand[0] = item->strand[0];
	pmat->array[i] = pbw;
	pmat->matrix[i] = pbw->data;
    }
    return pmat;
}

void perBaseMatrixAddOrigRegions(struct perBaseMatrix *pbm, struct bed6 *orig_regions)
/* add the original bed information for retrieval later.  it's really important */
/* that this is in the same order as the pbws in the array  */
{
    struct bed6 *region;
    int i;
    int num_regions = slCount(orig_regions);
    if (pbm->nrow != num_regions)
	errAbort("something's wrong");
    for (i = 0, region = orig_regions; ((i < pbm->nrow) && (region != NULL)); i++, region = region->next)
	pbm->array[i]->orig_bed = region;
}

void free_perBaseMatrix(struct perBaseMatrix **ppMat)
/* freeing up a wig matrix */
{
    struct perBaseMatrix *pmat = *ppMat;
    int i;
    for (i = 0; i < pmat->nrow; i++)
    {
	struct perBaseWig *pbw = pmat->array[i];
	perBaseWigFree(&pbw);
    }
    freeMem(pmat->array);
    freeMem(pmat->matrix);
    freez(ppMat);
}

static void wigOutJustVal(FILE *output, int decimals, double val)
{
    int dec = 3;
    if (decimals >= 0)
	dec = decimals;
    if (isnan(val))
	fprintf(output, "NA\n");
    else
	fprintf(output, "%0.*f\n", dec, val);
}

static void perBaseWigOutSection(struct perBaseWig *pbw, FILE *output, enum wigOutType wot,
				 int start, int end, int decimals, boolean nohead, boolean condense)
/* Simply output wigs in different ways.  */
{
    int i;
    if (wot == bedGraphOut)
    {
	if (condense)
	{
	    int j;
	    i = start;
	    while (i < end)
	    {
		j = i+1;
		while ((j < end) && (pbw->data[i] == pbw->data[j]))
		    j++;
		fprintf(output, "%s\t%d\t%d\t", pbw->chrom, pbw->chromStart + i, pbw->chromStart + j);
		wigOutJustVal(output, decimals, pbw->data[i]);
		i = j;
	    }
	}
	else
	{
	    for (i = start; i < end; i++)
	    {
		fprintf(output, "%s\t%d\t%d\t", pbw->chrom, pbw->chromStart + i, pbw->chromStart + i + 1);
		wigOutJustVal(output, decimals, pbw->data[i]);
	    }
	}
    }	    
    else if (wot == varStepOut)
    {
	if (!nohead)
	    fprintf(output, "variableStep chrom=%s span=1\n", pbw->chrom);
	for (i = start; i < end; i++)
	{
	    fprintf(output, "%d\t", pbw->chromStart + i + 1);
	    wigOutJustVal(output, decimals, pbw->data[i]);
	}
    }
    else if (wot == fixStepOut)
    {
	if (!nohead)
	    fprintf(output, "fixedStep chrom=%s start=%d step=1 span=1\n", pbw->chrom, pbw->chromStart + start + 1);
	for (i = start; i < end; i++)
	    wigOutJustVal(output, decimals, pbw->data[i]);
    }
    else
	errAbort("internal error: bad type of wig specified for output");
}
		

void perBaseWigOutput(struct perBaseWig *pbwList, FILE *output, enum wigOutType wot, 
		      int decimals, char *track_line, boolean nohead, boolean condense)
/* Output wiggle. */ 
{
    struct perBaseWig *pbw;
    if (track_line && !nohead)
	fprintf(output, "%s\n", track_line);
    for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
    {
	if (pbw->subsections)
	{
	    struct bed *bed;
	    for (bed = pbw->subsections; bed != NULL; bed = bed->next)
		perBaseWigOutSection(pbw, output, wot, bed->chromStart, bed->chromEnd, decimals, nohead, condense); 
	}
	else
	    perBaseWigOutSection(pbw, output, wot, 0, pbw->chromEnd-pbw->chromStart, decimals, nohead, condense);
    }
}

static void perBaseWigOutSectionNASkip(struct perBaseWig *pbw, FILE *output, enum wigOutType wot, int decimals, int start, int end,
				       boolean nohead, boolean condense)
/* Repeatedly outputs a section while skipping over NA data. */
{
    int i = start;
    while (i < end)
    {
	int s = 0;
	int e = 0;
	/* skip NA data */
	while ((i < end) && (isnan(pbw->data[i])))
	    i++;
	s = i;
	/* keep actual data */
	while ((i < end) && (!isnan(pbw->data[i])))
	    i++;
	e = i;
	if (s < end)
	    perBaseWigOutSection(pbw, output, wot, s, e, decimals, nohead, condense);
    }
}

void perBaseWigOutputNASkip(struct perBaseWig *pbwList, FILE *output, enum wigOutType wot, 
			    int decimals, char *track_line, boolean nohead, boolean condense)
/* Output wiggle skipping a specific value designated as NA. */ 
{
    struct perBaseWig *pbw;
    if (track_line && !nohead)
	fprintf(output, "%s\n", track_line);
    for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
    {
	if (pbw->subsections)
	{
	    struct bed *bed;
	    for (bed = pbw->subsections; bed != NULL; bed = bed->next)
		perBaseWigOutSectionNASkip(pbw, output, wot, decimals, bed->chromStart, 
					   bed->chromEnd, nohead, condense); 
	}
	else
	    perBaseWigOutSectionNASkip(pbw, output, wot, decimals, 0, pbw->chromEnd-pbw->chromStart, nohead, condense);
    }
}

static int pos_compare(const void *a, const void *b)
/* used by qsort for the starts */
{
    return *(int *)a - *(int *)b;
}

static int pos_compare2(const void *a, const void *b)
/* used by qsort for the middles. sorts tuples */
{
    int *pA = (int *)a;
    int *pB = (int *)b;
    return pA[0] - pB[0];
}

static int unique_nums(const int *array, int size)
{
    int i;
    int unique = 0;
    int last = -1;
    for (i = 0; i < size; i++)
    {
	if (array[i] != last)
	    unique++;
	last = array[i];
    }
    return unique;
}

static int unique_nums2(const int *array, int size)
{
    int i;
    int unique = 0;
    int last = -1;
    for (i = 0; i < size * 2; i += 2)
    {
	if (array[i] != last)
	    unique++;
	last = array[i];
    }
    return unique;
}

static void fold_array(int *big_array, int big_array_size, int *unique_array, int *counts_array)
/* convert array of starts into unique array of starts + array of counts */
{
    int i = 0; 
    int j = -1;
    while (i < big_array_size)
    {
	j++;
	unique_array[j] = big_array[i];
	counts_array[j] = 0;
	while ((i < big_array_size) && (big_array[i] == unique_array[j]))
	{
	    i++;
	    counts_array[j]++;
	}
    }
}

static void fold_array2(int *big_array, int big_array_size, int *unique_array, int *counts_array, int *ave_sizes)
/* convert array of starts into unique array of starts + array of counts */
{
    int i = 0; 
    int j = -1;
    while (i < big_array_size)
    {
	j++;
	unique_array[j] = big_array[i];
	counts_array[j] = 0;
	while ((i < big_array_size) && (big_array[i] == unique_array[j]))
	{
	    counts_array[j]++;
	    ave_sizes[j] += big_array[i+1];
	    i += 2;
	}
	ave_sizes[j] /= counts_array[j];
    }
}

struct middles *beds_to_middles(struct bed6 *bedList)
/* minimal information for stacking middles */
{
    struct bed6 *cur;
    struct middles *mids;
    int *pos_mids;
    int size = slCount(bedList);
    int i = 0;
    if (size == 0)
	return NULL;
    AllocArray(pos_mids, size * 2);
    for (cur = bedList; (i < size) && (cur != NULL); i++, cur = cur->next)
    {	
	pos_mids[i*2] = (cur->chromStart + cur->chromEnd)/2;
	pos_mids[i*2 + 1] = cur->chromEnd - cur->chromStart;
	/* note: can be out-of-range positions. we don't deal with this here */
    }
    qsort(pos_mids, size, sizeof(int) * 2, pos_compare2);
    AllocVar(mids);
    mids->num_mids = unique_nums2(pos_mids, size);
    AllocArray(mids->mids, mids->num_mids);
    AllocArray(mids->counts, mids->num_mids);
    AllocArray(mids->ave_sizes, mids->num_mids);
    fold_array2(pos_mids, size * 2, mids->mids, mids->counts, mids->ave_sizes);
    freeMem(pos_mids);
    return mids;
}

static struct starts *beds_to_starts(struct bed6 *bedList, boolean both_ends)
/* convert bigBedIntervals to minimal information needed for phasing */
{
    struct bed6 *cur;
    struct starts *starts;
    int *pos_starts = NULL;
    int *neg_starts = NULL;
    int total_size = slCount(bedList);
    int num_ps = 0;
    int num_ns = 0;
    /* first build arrays to hold all the starts and make sure it's sorted */
    /* on the minus strand */
    if (total_size == 0)
	return NULL;
    AllocArray(pos_starts, total_size*2);
    AllocArray(neg_starts, total_size);
    for (cur = bedList; cur != NULL; cur = cur->next)
    {
	if (both_ends)
	{
	    pos_starts[num_ps++] = cur->chromStart;
	    neg_starts[num_ns++] = cur->chromEnd;	    
	}
	else if (cur->strand[0] == '+')
	    pos_starts[num_ps++] = cur->chromStart;
	else
	    neg_starts[num_ns++] = cur->chromEnd;
    }
    if (both_ends)
	qsort(pos_starts, num_ps, sizeof(int), pos_compare);
    else
	qsort(neg_starts, num_ns, sizeof(int), pos_compare);
    /* find the number of unique starts */
    AllocVar(starts);
    starts->num_pos_starts = unique_nums(pos_starts, num_ps);
    starts->num_neg_starts = unique_nums(neg_starts, num_ns);
    /* init arrays */
    if (starts->num_pos_starts > 0)
    {
	AllocArray(starts->pos_starts, starts->num_pos_starts);
	AllocArray(starts->pos_counts, starts->num_pos_starts);
    }
    if (starts->num_neg_starts > 0)
    {
	AllocArray(starts->neg_starts, starts->num_neg_starts);
	AllocArray(starts->neg_counts, starts->num_neg_starts);
    }
    /* convert the starts into unique starts with counts */
    fold_array(pos_starts, num_ps, starts->pos_starts, starts->pos_counts);
    fold_array(neg_starts, num_ns, starts->neg_starts, starts->neg_counts);
    freeMem(pos_starts);
    freeMem(neg_starts);
    return starts;
}

struct starts *metaBig_get_starts(struct metaBig *mb, char *chrom, unsigned start, unsigned end)
/* return starts struct used with the phasogram/distogram programs */
{
    struct starts *starts = NULL;
    struct lm *lm = lmInit(0);
    struct bed6 *bedList = metaBigBed6Fetch(mb, chrom, start, end, lm);
    starts = beds_to_starts(bedList, FALSE);
    lmCleanup(&lm);
    return starts;
}

struct starts *metaBig_get_starts_both_ends(struct metaBig *mb, char *chrom, unsigned start, unsigned end)
/* return starts struct used with the phasogram/distogram programs */
{
    struct starts *starts = NULL;
    struct lm *lm = lmInit(0);
    struct bed6 *bedList = metaBigBed6Fetch(mb, chrom, start, end, lm);
    starts = beds_to_starts(bedList, TRUE);
    lmCleanup(&lm);
    return starts;
}

struct middles *metaBig_get_middles(struct metaBig *mb, char *chrom, unsigned start, unsigned end)
/* return middles struct used with stackreads program */
{
    struct middles *mids = NULL;
    struct lm *lm = lmInit(0);
    struct bed6 *bedList = metaBigBed6Fetch(mb, chrom, start, end, lm);
    mids = beds_to_middles(bedList);
    lmCleanup(&lm);
    return mids;
}

static boolean local_file(char *filename)
/* return TRUE if the file is avialable locally */
{
    char *s = cloneString(filename);
    char *colon = strchr(s, ':');
    boolean ret = FALSE;
    if (colon)
	*colon = '\0';
    if (fileExists(s))
	ret = TRUE;
    freeMem(s);
    return ret;
}

struct bed *metaBig_chopGenome(struct metaBig *mb, int size)
/* return a bed of regularly-sized intervals (given) from the chromSizeHash */
{
    struct hashEl *list = hashElListHash(mb->chromSizeHash);
    struct hashEl *el;
    struct bed *bed_list = NULL;
    for (el = list; el != NULL; el = el->next)
    {
	char *chrom = el->name;
	int chromSize = ptToInt(el->val);
	int i;
	for (i = 0; i < chromSize; i += size)
	{
	    struct bed *newbed;
	    AllocVar(newbed);
	    newbed->chrom = cloneString(chrom);
	    newbed->chromStart = i;
	    if ((i + size) > chromSize)
		newbed->chromEnd = chromSize;
	    else
		newbed->chromEnd = i+size;
	    newbed->name = cloneString(".");
	    slAddHead(&bed_list, newbed);
	}
    }
    slReverse(&bed_list);
    return bed_list;
}

static void startsFree(struct starts **pStarts)
/* free a starts */
{
    struct starts *starts = *pStarts;
    freeMem(starts->pos_starts);
    freeMem(starts->neg_starts);
    freez(pStarts);
}

void startsFreeList(struct starts **pList)
/* free a list of starts */
{
    struct starts *reg;
    while ((reg = slPopHead(pList)) != NULL)
	startsFree(&reg);
    freez(pList);
}

void middlesFreeList(struct middles **pList)
/* free a list of middles */
{
    struct middles *mid;
    while ((mid = slPopHead(pList)) != NULL)
    {
	freeMem(mid->mids);
	freeMem(mid->counts);
	freez(&mid);
    }
    freez(pList);
}
