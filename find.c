/* bwtool_find - find parts of the  data to  */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <jkweb/common.h>
#include <jkweb/linefile.h>
#include <jkweb/hash.h>
#include <jkweb/options.h>
#include <jkweb/sqlNum.h>
#include <jkweb/basicBed.h>
#include <jkweb/bigWig.h>
#include <beato/bigs.h>
#include "bwtool.h"
#include "bwtool_shared.h"
#include <jkweb/rangeTree.h>
#include <beato/extrema.h>

#include <float.h>

void usage_find()
/* Explain usage and exit. */
{
errAbort(
  "bwtool find - find parts of genome that satisfy given constraints\n"
  "   in the bigWig data.\n"
  "usage:\n"
  "   bwtool find <operator> [operator options] input.bw[:chr:start-end] output.bed\n"
  "operators:\n"
  "   local-extrema          output is 6-field bedGraph with +/- indicating if the extrema\n"
  "                          is a local minimum or maximum.\n"
  "        -min-sep=d        to avoid local clusters of extrema, define a minimum distance\n"
  "                          they must be separated by (d)\n"
  "        -maxima           only find local maxima\n"
  "        -minima           only find local minima\n"
  "   less|less-equal|more|more-equal|equal|not <number>\n"
  "   maxima <regions.bed>   find the highest value in each region given in a bed and output\n"
  "                          the same bed as bed12.\n"
  "        -with-max         output the maximum value just after the bed12.\n"
  "        -median-base      if there are multiple tied maxima, output the median base location\n"
  "                          instead of all of them.\n"
  );
}

static enum ex_removal get_removal(struct hash *options)
/* convert -minima and -maxima options to enum */
{
    boolean maxima = (hashFindVal(options, "maxima") != NULL) ? TRUE : FALSE;
    boolean minima = (hashFindVal(options, "minima") != NULL) ? TRUE : FALSE;
    if (maxima && minima)
	errAbort("cannot specify both -maxima and -minima");
    if (maxima)
	return remove_min;
    if (minima)
	return remove_max;
    return no_removal;
}

void bwtool_find_extrema(struct hash *options, char *favorites, char *regions, unsigned decimals, double fill, char *bigfile, char *tmp_dir, char *outputfile)
/* find local extrema */
{
    unsigned min_sep = sqlUnsigned((char *)hashOptionalVal(options, "min-sep", "0"));
    char *other_bigfile = (char *)hashOptionalVal(options, "against", NULL);
    enum ex_removal rem = get_removal(options);
    struct metaBig *main_big = metaBigOpen_check(bigfile, tmp_dir, regions);
    struct metaBig *other_big = NULL;
    struct extrema *main_list;
    struct extrema *other_list = NULL;
    struct extrema *ex;
    unsigned shift = 0;
    FILE *out;
    if (other_bigfile)
    {
	char *num;
	if (rem == no_removal)
	    errAbort("must specify either -maxima or -minima with -against");
	if (!strchr(other_bigfile, ','))
	    errAbort("must specify shift limit in -against option");
	num = chopPrefixAt(other_bigfile, ',');
	shift = sqlUnsigned(num);
	other_big = metaBigOpen_check(other_bigfile, tmp_dir, regions);
    }
    if (!main_big || (!other_big && other_bigfile))
	errAbort("could not open bigWig file");
    main_list = extrema_find(main_big, min_sep, rem);
    slReverse(&main_list);
    if (other_bigfile)
    {
	other_list = extrema_find(other_big, min_sep, rem);
	extrema_find_shifts(main_list, other_list, shift);
    }
    metaBigClose(&main_big);
    if (other_big)
	metaBigClose(&other_big);
    out = mustOpen(outputfile, "w");
    if (other_bigfile)
	for (ex = main_list; ex != NULL; ex = ex->next)
	    fprintf(out, "%s\t%d\t%d\t%d\t1000\t%c\n", ex->chrom, ex->chromStart, ex->chromStart+1, (int)ex->val, ex->min_or_max);
    else
    {
	slSort(&main_list, extrema_bed_cmp);
	for (ex = main_list; ex != NULL; ex = ex->next)
	    fprintf(out, "%s\t%d\t%d\t%0.*f\t1000\t%c\n", ex->chrom, ex->chromStart, ex->chromStart+1, decimals, ex->val, ex->min_or_max);
    }
    carefulClose(&out);
    extrema_free_list(&main_list);
}

static boolean fit_thresh(double val, double thresh, enum bw_op_type op)
/* just a big switch depending on the operator */
{
    boolean ret = FALSE;
    switch (op)
    {
    case less:
    {
	if (val < thresh)
	    ret = TRUE;
	break;
    }
    case less_equal:
    {
	if (val <= thresh)
	    ret = TRUE;
	break;
    }
    case more:
    {
	if (val > thresh)
	    ret = TRUE;
	break;
    }
    case more_equal:
    {
	if (val >= thresh)
	    ret = TRUE;
	break;
    }
    case equal:
    {
	if (val == thresh)
	    ret = TRUE;
	break;
    }
    case not_equal:
    {
	if (val != thresh)
	    ret = TRUE;
	break;
    }
    case invalid:
    default:
    {
	ret = FALSE;
    }}
    return ret;
}

void bwtool_find_thresh(struct hash *options, char *favorites, char *regions, double fill,
			char *thresh_type, char *thresh_s, char *bigfile, char *tmp_dir, char *outputfile)
/* the other kind of finding, based on thresholding. */
{
    boolean inverse = (hashFindVal(options, "inverse") != NULL) ? TRUE : FALSE;
    enum bw_op_type op= get_bw_op_type(thresh_type, inverse);
    struct metaBig *mb = metaBigOpen_check(bigfile, tmp_dir, regions);
    double thresh = sqlDouble(thresh_s);
    FILE *out = mustOpen(outputfile, "w");
    struct bed out_bed;
    struct bed *section;
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbwList = perBaseWigLoadContinue(mb, section->chrom, section->chromStart,
							      section->chromEnd);
	struct perBaseWig *pbw;
	int i, len;
	if (pbwList)
	{
	    out_bed.chrom = pbwList->chrom;
	    for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
	    {
		i = 0;
		len = pbw->chromEnd - pbw->chromStart;
		out_bed.chromStart = out_bed.chromEnd = 0;
		while (i < len)
		{
		    while ((i < len) && (!fit_thresh(pbw->data[i], thresh, op)))
			i++;
		    out_bed.chromStart = i + pbw->chromStart;
		    while ((i < len) && (fit_thresh(pbw->data[i], thresh, op)))
			i++;
		    out_bed.chromEnd = i + pbw->chromStart;
		    if (out_bed.chromEnd > out_bed.chromStart)
			bedTabOutN(&out_bed, 3, out);
		}
	    }
	perBaseWigFree(&pbwList);
	}
    }
    metaBigClose(&mb);
    carefulClose(&out);
}

static struct bed *bed12FromBed6(struct bed6 **pList)
/* Destroy bed6 list in favor of a typical bed */
{
    struct bed *list = NULL;
    struct bed6 *top;
    while ((top = slPopHead(pList)) != NULL)
    {
	struct bed *uno;
	AllocVar(uno);
	uno->chrom = cloneString(top->chrom);
	uno->chromStart = top->chromStart;
	uno->chromEnd = top->chromEnd;
	uno->name = cloneString(top->name);
	uno->score = top->score;
	uno->strand[0] = top->strand[0];
	uno->strand[1] = '\0';
	uno->thickStart = uno->chromStart;
	uno->thickEnd = uno->chromEnd;
	slAddHead(&list, uno);
	bed6Free(&top);
    }
    freez(pList);
    slReverse(&list);
    return list;
}

static int median_base_calc(struct slInt **pList)
/* in a list of ints, find the median.  straight-forward enough. */
{
    struct slInt *list;
    if (!pList || !*pList)
	return -1;
    int size = slCount(*pList);
    list = *pList;
    if (size == 1)
	return list->val;
    int i = 0;
    struct slInt *cur;
    slSort(pList, slIntCmp);
    list = *pList;
    if (size % 2 == 0)
    {	
	for (cur = list, i = 0; ((cur != NULL) && (i < (size/2)-1)); cur = cur->next, i++)
	     ;
	return (cur->val + cur->next->val)/2;
    } 
    for (cur = list, i = 0; ((cur != NULL) && (i < (size/2))); cur = cur->next, i++)
	;
    return cur->val;
}

void bwtool_find_max(struct hash *options, char *favorites, char *regions, double fill,
		     char *bigfile, char *tmp_dir, char *outputfile)
/* find max points in a range */
{
    boolean med_base = (hashFindVal(options, "median-base") != NULL) ? TRUE : FALSE;
    boolean with_max = (hashFindVal(options, "with-max") != NULL) ? TRUE : FALSE;
    struct metaBig *mb = metaBigOpen_check(bigfile, tmp_dir, NULL);
    FILE *out = mustOpen(outputfile, "w");
    struct bed6 *sections6 = readBed6Soft(regions);
    struct bed *sections = bed12FromBed6(&sections6);
    struct bed *section;
    for (section = sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbwList = perBaseWigLoadContinue(mb, section->chrom, section->chromStart,
							      section->chromEnd);
	struct perBaseWig *pbw;
	struct slInt *ii;
	int i, size;
	double max = -DBL_MAX;
	struct slInt *list = NULL;
	for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
	{
	    int pbw_off = pbw->chromStart - section->chromStart;
	    for (i = 0; i < pbw->len; i++)
	    {
		if (pbw->data[i] > max)
		{
		    slFreeList(&list);
		    struct slInt *new_int = slIntNew(i + pbw_off);
		    slAddHead(&list, new_int);
		    max = pbw->data[i];
		}
		else if (pbw->data[i] == max)
		{
		    struct slInt *new_int = slIntNew(i + pbw_off);
		    slAddHead(&list, new_int);
		}
	    }
	}
	slReverse(&list);
	if (list)
	{
	    size = slCount(list);
	    if (med_base)
	    {
		section->blockCount = 1;
		AllocArray(section->blockSizes, sizeof(int));
		AllocArray(section->chromStarts, sizeof(int));
		section->blockSizes[0] = 1;
		section->chromStarts[0] = median_base_calc(&list);
	    }
	    else
	    {
		section->blockCount = size;
		AllocArray(section->blockSizes, sizeof(int) * size);
		AllocArray(section->chromStarts, sizeof(int) * size);
		for (i = 0, ii = list; (i < size) && (ii != NULL); i++, ii = ii->next)
		{
		    section->blockSizes[i] = 1;
		    section->chromStarts[i] = ii->val;
		}
	    }
	    if (!with_max)
		bedTabOutN(section, 12, out);
	    else
	    {
		bedOutputN(section, 12, out, '\t', '\t');
		fprintf(out, "%f\n", max);
	    }
	    slFreeList(&list);
	}
	perBaseWigFree(&pbwList);
    }
    metaBigClose(&mb);
    bedFreeList(&sections);
    carefulClose(&out);
}
