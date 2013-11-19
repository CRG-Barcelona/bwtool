/* bwtool_find - find parts of the  data to  */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif 

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "basicBed.h"
#include "bigWig.h"
#include "bigs.h"
#include "bwtool.h"
#include "bwtool_shared.h"
#include "rangeTree.h"
#include "extrema.h"

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
  "   local-extrema        output is 6-field bedGraph with +/- indicating if the extrema\n"
  "                        is a local minimum or maximum.\n"
  "        -min-sep=d      to avoid local clusters of extrema, define a minimum distance\n"
  "                        they must be separated by (d)\n"
  "        -maxima         only find local maxima\n"
  "        -minima         only find local minima\n"
  "   less|less-equal|more|more-equal|equal|not <number>\n"
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

void bwtool_find_extrema(struct hash *options, char *favorites, char *regions, unsigned decimals, double fill, char *bigfile, char *outputfile)
/* find local extrema */
{
    unsigned min_sep = sqlUnsigned((char *)hashOptionalVal(options, "min-sep", "0"));
    char *other_bigfile = (char *)hashOptionalVal(options, "against", NULL);
    enum ex_removal rem = get_removal(options);
    struct metaBig *main_big = metaBigOpen_check(bigfile, regions);
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
	other_big = metaBigOpen_check(other_bigfile, regions);
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
	for (ex = main_list; ex != NULL; ex = ex->next)
	    fprintf(out, "%s\t%d\t%d\t%0.*f\t1000\t%c\n", ex->chrom, ex->chromStart, ex->chromStart+1, decimals, ex->val, ex->min_or_max);
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
			char *thresh_type, char *thresh_s, char *bigfile, char *outputfile)
/* the other kind of finding, based on thresholding. */
{
    boolean inverse = (hashFindVal(options, "inverse") != NULL) ? TRUE : FALSE;
    enum bw_op_type op= get_bw_op_type(thresh_type, inverse);
    struct metaBig *mb = metaBigOpen_check(bigfile, regions);
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
    metaBigClose(&mb);
    carefulClose(&out);
}
