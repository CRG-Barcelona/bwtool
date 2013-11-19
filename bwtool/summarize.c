/* bwtool_chromgraph - chromgraph plotting  */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif 

#include "common.h"
#include "obscure.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "basicBed.h"
#include "bigWig.h"
#include "bigs.h"
#include "stuff.h"
#include "bwtool.h"
#include "bwtool_shared.h"

#define NANUM sqrt(-1)

void usage_summary()
/* Explain usage of the summarize program and exit. */
{
errAbort(
  "bwtool summary - provide some summary stats for each region in a bed file\n"
  "   or at regular intervals.\n"
  "usage:\n"
  "   bwtool summary loci input.bw[:chr:start-end] output.txt\n" 
  "where:\n"
  "   -\"loci\" corresponds to either (a) a bed file with regions to summarize or\n"
  "    (b) a size of interval to summarize genome-wide.\n"
  "options:\n"
  "   -with-quantiles  output 10%%/25%%/75%%/90%% quantiles as well surrounding the\n"
  "                    median.  With -total, this essentially provides a boxplot.\n"
  "   -with-sum-of-squares\n"
  "                    output sum of squared deviations from the mean along with \n"
  "                    the other fields\n"
  "   -with-sum        output sum, also\n"
  "   -keep-bed        if the loci bed is given, keep as many bed file\n"
  "   -total           only output a summary as if all of the regions are pasted\n"
  "                    together\n"
  "   -header          put in a header (fields are easy to forget)\n"
  );
}

void fill_na(double *data, int size, unsigned decimals)
{
    double na = NANUM;
    int i = 0;
    for (i = 0; i < size; i++)
    {
	double epsilon = 1/pow(10, decimals+1);
	if ((data[i] < epsilon) && (data[i] > -1*epsilon))
	    data[i] = na;
    }
}

double sumOfSquares(int num_data, double *vector, double mean)
/* return sum of squares */
{
    int i;
    double ret = 0;
    for (i = 0; i < num_data; i++)
	ret += pow(vector[i] - mean, 2);
    return ret;
}

void summary_loop(struct perBaseWig *pbw, unsigned decimals, FILE *out, struct bed *section, int bed_size, boolean use_rgb, boolean zero_remove, boolean with_quants, boolean with_sos, boolean with_sum)
/* at each iteration of */
{
    int i;
    int size = pbw->chromEnd - pbw->chromStart;
    double *vector = pbw->data;
    if (zero_remove)
	fill_na(vector, size, decimals);
    unsigned num_data = doubleWithNASort(size, vector);
    if (num_data > 0)
    {
	double mean = 0;
	double sum = 0;
	double min = DBL_MAX;
	double max = -1 * DBL_MAX;
	double median = doubleWithNAMedianAlreadySorted(num_data, vector);
	double first_10p = -1;
	double first_quart = -1;
	double third_quart = -1;
	double last_10p = -1;
	double sos = -1;
	if (with_quants)
	{
	    first_10p = doubleWithNAInvQuantAlreadySorted(num_data, vector, 10, TRUE);
	    first_quart = doubleWithNAInvQuantAlreadySorted(num_data, vector, 4, TRUE);
	    third_quart = doubleWithNAInvQuantAlreadySorted(num_data, vector, 4, FALSE);
	    last_10p = doubleWithNAInvQuantAlreadySorted(num_data, vector, 10, FALSE);
	}
	for (i = 0; i < num_data; i++)
	{
	    sum += vector[i];
	    if (vector[i] > max)
		max = vector[i];
	    if (vector[i] < min)
		min = vector[i];
	}
	mean = sum/num_data;
	if (with_sos)
	    sos = sumOfSquares(num_data, vector, mean);
	bedOutFlexible(section, bed_size, out, '\t', '\t', use_rgb);
	if (with_quants)
	    fprintf(out, "%d\t%d\t%0.*f\t%0.*f\t%0.*f\t%0.*f\t%0.*f\t%0.*f\t%0.*f\t%0.*f", size, 
		    num_data, decimals, min, decimals, max, decimals, mean, decimals, first_10p, decimals, first_quart, decimals, median, decimals, third_quart, decimals, last_10p);
	else
	    fprintf(out, "%d\t%d\t%0.*f\t%0.*f\t%0.*f\t%0.*f", size, 
		num_data, decimals, min, decimals, max, decimals, mean, decimals, median);
	if (with_sos)
	    fprintf(out, "\t%0.*f", decimals, sos);
	if (with_sum)
	    fprintf(out, "\t%0.*f", decimals, sum);
	fprintf(out, "\n");
    }
    else
    {
	bedOutFlexible(section, bed_size, out, '\t', '\t', use_rgb);
	if (with_quants)
	    fprintf(out, "%d\tna\tna\tna\tna\tna\tna\tna\tna\tna", size);
	else
	    fprintf(out, "%d\tna\tna\tna\tna\tna", size);
	if (with_sos)
	    fprintf(out, "\tna");
	if (with_sum)
	    fprintf(out, "\tna");
	fprintf(out, "\n");
    }
}

void bwtool_summary_bed(struct metaBig *mb, unsigned decimals, struct bed *bed_list, int bed_size, boolean use_rgb, FILE *out, double fill, boolean zero_remove, boolean with_quants, boolean with_sos, boolean with_sum, boolean total)
/* if the "loci" ends up being a bed file */
{
    if (total)
    {
	struct perBaseWig *big_pbw = perBaseWigLoadHuge(mb, bed_list);
	struct bed *big_bed;	
	AllocVar(big_bed);
	big_bed->chrom = cloneString(big_pbw->chrom);
	big_bed->chromStart = big_pbw->chromStart;
	big_bed->chromEnd = big_pbw->chromEnd;
	summary_loop(big_pbw, decimals, out, big_bed, 3, FALSE, zero_remove, with_quants, with_sos, with_sum);
	bedFree(&big_bed);
    }
    else
    {
	struct bed *section;
	for (section = bed_list; section != NULL; section = section->next)
	{
	    struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, section->chrom, section->chromStart, section->chromEnd, FALSE, fill);
	    summary_loop(pbw, decimals, out, section, bed_size, use_rgb, zero_remove, with_quants, with_sos, with_sum);
	    perBaseWigFreeList(&pbw);
	}
    }
}

void bwtool_summary(struct hash *options, char *favorites, char *regions, unsigned decimals,
		    double fill, char *loci_s, char *bigfile, char *outputfile)
/* bwtool_summary - main for the summarize program */
{
    boolean zero_remove = (hashFindVal(options, "zero-remove") != NULL) ? TRUE : FALSE;
    boolean keep_bed = (hashFindVal(options, "keep-bed") != NULL) ? TRUE : FALSE;
    boolean with_quants = (hashFindVal(options, "with-quantiles") != NULL) ? TRUE : FALSE;
    boolean with_sos = (hashFindVal(options, "with-sum-of-squares") != NULL) ? TRUE : FALSE;
    boolean with_sum = (hashFindVal(options, "with-sum") != NULL) ? TRUE : FALSE;
    boolean total = (hashFindVal(options, "total") != NULL) ? TRUE : FALSE;
    struct metaBig *mb = metaBigOpen_check(bigfile, regions);
    double filll = (hashFindVal(options, "zero-fill") != NULL) ? 0 : fill;
    if ((fill == 0) && total)
	errAbort("-total incompatible with -zero-fill");
    else if (!isnan(filll) && total)
	errAbort("-total incompatible with -fill");
    if (total && keep_bed)
	warn("-keep-bed useless with -total");
    if (mb->type != isaBigWig)
	errAbort("file not bigWig type");
    FILE *out = mustOpen(outputfile, "w");
    boolean header = (hashFindVal(options, "header") != NULL) ? TRUE : FALSE;
    struct bed *bed_list = NULL;
    boolean use_rgb = FALSE;
    int bed_size = 3;
    if (fileExists(loci_s))
	bedLoadAllReturnFieldCountAndRgbAtLeast3(loci_s, &bed_list, &bed_size, &use_rgb);
    else
    {
	unsigned interval = sqlUnsigned(loci_s);
	bed_list = metaBig_chopGenome(mb, interval);
    }
    if (!keep_bed)
	bed_size = 3;
    if (header)
    {
	int i;
	fprintf(out, "#chrom\tstart\tend");
	if (bed_size > 3)
	    fprintf(out, "\tname");
	if (bed_size > 4)
	    fprintf(out, "\tscore");
	if (bed_size > 5)
	    fprintf(out, "\tstrand");
	if (bed_size > 6)
	    fprintf(out, "\tthick_start");
	if (bed_size > 7)
	    fprintf(out, "\tthick_end");
	if ((bed_size > 8) && (use_rgb))
	    fprintf(out, "\trgb");
	if ((bed_size > 8) && (!use_rgb))
	    fprintf(out, "\treserved");
	if (bed_size > 9)
	    fprintf(out, "\tblocks");
	if (bed_size > 10)
	    fprintf(out, "\tblock_sizes");
	if (bed_size > 11)
	    fprintf(out, "\tblock_starts");
	for (i = 12; i < bed_size; i++)
	    fprintf(out, "\tunknown_field");
	if (with_quants)
	    fprintf(out, "\tsize\tnum_data\tmin\tmax\tmean\tfirst_10%%\t1st_quart\tmedian\t3rd_quart\tlast_10%%");
	else
	    fprintf(out, "\tsize\tnum_data\tmin\tmax\tmean\tmedian");
	if (with_sos)
	    fprintf(out, "\tsum_of_squares");
	if (with_sum)
	    fprintf(out, "\tsum");
	fprintf(out, "\n");
    }
    bwtool_summary_bed(mb, decimals, bed_list, bed_size, use_rgb, out, fill, zero_remove, with_quants, with_sos, with_sum, total);
    bedFreeList(&bed_list);
    carefulClose(&out);
    metaBigClose(&mb);
}
