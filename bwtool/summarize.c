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

void summary_loop(struct metaBig *mb, unsigned decimals, FILE *out, char *chrom, unsigned chromStart, unsigned chromEnd, boolean zero_fill)
/* at each iteration of */
{
    int i;
    struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, chrom, chromStart, chromEnd, FALSE);
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
	fprintf(out, "%s\t%d\t%d\t%d\t%d\t%0.*f\t%0.*f\t%0.*f\t%0.*f\n", chrom, chromStart, chromEnd, size, 
		num_data, decimals, min, decimals, max, decimals, mean, decimals, median);
    }
    else
    {
	fprintf(out, "%s\t%d\t%d\t0\tna\tna\tna\tna\tna\n", chrom, chromStart, chromEnd, size);
    }
    perBaseWigFreeList(&pbw);
}

void bwtool_summary_bed(struct metaBig *mb, unsigned decimals, char *bedfile, FILE *out, boolean zero_fill)
/* if the "loci" ends up being a bed file */
{
    struct bed *section;
    struct bed *summary_regions = bedLoadNAll(bedfile, 3);
    for (section = summary_regions; section != NULL; section = section->next)
	summary_loop(mb, decimals, out, section->chrom, section->chromStart, section->chromEnd, zero_fill);
    bedFreeList(&summary_regions);
}

void bwtool_summary_interval(struct metaBig *mb, unsigned decimals, unsigned interval, FILE *out, boolean zero_fill)
/* if the "loci" ends up being an interval */
{
    struct hashEl *chroms = hashElListHash(mb->chromSizeHash);
    struct hashEl *el;
    for (el = chroms; el != NULL; el = el->next)
    {
	int size = ptToInt(el->val);
	unsigned start = 0;
	unsigned end;
	for (start = 0; start < size; start += interval)
	{
	    end = start + interval;
	    if (end > size)
		end = size;
	    summary_loop(mb, decimals, out, el->name, start, end, zero_fill);
	}
    }
    hashElFreeList(&chroms);
}

void bwtool_summary(struct hash *options, char *favorites, char *regions, unsigned decimals,
		    char *loci_s, char *bigfile, char *outputfile)
/* bwtool_summary - main for the summarize program */
{
    boolean zero_fill = (hashFindVal(options, "zero-fill") != NULL) ? TRUE : FALSE;
    struct metaBig *mb = metaBigOpen_favs(bigfile, regions, favorites);
    if (mb->type != isaBigWig)
	errAbort("file not bigWig type");
    FILE *out = mustOpen(outputfile, "w");
    boolean header = (hashFindVal(options, "header") != NULL) ? TRUE : FALSE;
    if (header)
	fprintf(out, "#chrom\tstart\tend\tsize\tnum_data\tmin\tmax\tmean\tmedian\n");
    if (fileExists(loci_s))
	bwtool_summary_bed(mb, decimals, loci_s, out, zero_fill);
    else
    {
	unsigned interval = sqlUnsigned(loci_s);
	bwtool_summary_interval(mb, decimals, interval, out, zero_fill);
    }
    carefulClose(&out);
    metaBigClose(&mb);
}
