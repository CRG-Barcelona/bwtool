/* bwtool_chromgraph - chromgraph plotting  */

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
  "   -zero-fill       treat NA regions as zero\n"
  "   -keep-bed        if the loci bed is given, keep as many bed file\n"
  "   -header          put in a header (fields are easy to forget)\n"
  );
}

void fill_zeros(double *data, int size)
/* with the -zero-fill option, add zeroes */
{
    int i;
    for (i = 0; i < size; i++)
	if (isnan(data[i]))
	    data[i] = 0;
}

void summary_loop(struct metaBig *mb, unsigned decimals, FILE *out, struct bed *section, int bed_size, boolean use_rgb, boolean zero_fill)
/* at each iteration of */
{
    int i;
    struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, section->chrom, section->chromStart, section->chromEnd, FALSE);
    int size = pbw->chromEnd - pbw->chromStart;
    double *vector = pbw->data;
    if (zero_fill)
	fill_zeros(vector, size);
    unsigned num_data = doubleWithNASort(size, vector);
    if (num_data > 0)
    {
	double mean = 0;
	double sum = 0;
	double min = DBL_MAX;
	double max = -1 * DBL_MAX;
	double median = doubleWithNAMedianAlreadySorted(num_data, vector);
	for (i = 0; i < num_data; i++)
	{
	    sum += vector[i];
	    if (vector[i] > max)
		max = vector[i];
	    if (vector[i] < min)
		min = vector[i];
	}
	mean = sum/num_data;
	bedOutFlexible(section, bed_size, out, '\t', '\t', use_rgb);
	fprintf(out, "%d\t%d\t%0.*f\t%0.*f\t%0.*f\t%0.*f\n", size, 
		num_data, decimals, min, decimals, max, decimals, mean, decimals, median);
    }
    else
    {
	bedOutFlexible(section, bed_size, out, '\t', '\t', use_rgb);
	fprintf(out, "0\tna\tna\tna\tna\tna\n");
    }
    perBaseWigFreeList(&pbw);
}

void bwtool_summary_bed(struct metaBig *mb, unsigned decimals, struct bed *bed_list, int bed_size, boolean use_rgb, FILE *out, boolean zero_fill)
/* if the "loci" ends up being a bed file */
{
    struct bed *section;
    for (section = bed_list; section != NULL; section = section->next)
	summary_loop(mb, decimals, out, section, bed_size, use_rgb, zero_fill);
}

void bwtool_summary(struct hash *options, char *favorites, char *regions, unsigned decimals,
		    char *loci_s, char *bigfile, char *outputfile)
/* bwtool_summary - main for the summarize program */
{
    boolean zero_fill = (hashFindVal(options, "zero-fill") != NULL) ? TRUE : FALSE;
    boolean keep_bed = (hashFindVal(options, "keep-bed") != NULL) ? TRUE : FALSE;
    struct metaBig *mb = metaBigOpen_favs(bigfile, regions, favorites);
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
	fprintf(out, "\tsize\tnum_data\tmin\tmax\tmean\tmedian\n");
    }
    bwtool_summary_bed(mb, decimals, bed_list, bed_size, use_rgb, out, zero_fill);
    bedFreeList(&bed_list);
    carefulClose(&out);
    metaBigClose(&mb);
}
