/* bwtool_distrib - distribution plotting  */

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "basicBed.h"
#include "bigWig.h"
#include "bigs.h"
#include "bwtool.h"

void usage_distrib()
/* Explain usage of distribution program and exit. */
{
errAbort(
  "bwtool distribution - produce plot data as the frequency of values\n"
  "   seen in the bigWig (converted to integers)\n"
  "usage:\n"
  "   bwtool distribution input.bw[:chr:start-end] output.txt\n" 
  "options:\n"
  "     -mult=m      multiply data by a number so the range is altered\n"
  );
}

void bwtool_distrib(struct hash *options, char *favorites, char *regions, unsigned decimals,
		    char *bigfile, char *outputfile)
/* bwtool_distrib - main for distribution program */
{
    struct metaBig *mb = metaBigOpen_favs(bigfile, regions, favorites);
    FILE *output = mustOpen(outputfile, "w");
    struct bed *section;
    int low = 0;
    int high = 0;
    int size = 0;
    int i, j;
    long *counts;
    struct bbiSummaryElement summary;
    if (mb->type != isaBigWig)
	errAbort("file not bigWig type");
    /* max - min will be the size of the count array */
    double m = sqlDouble((char *)hashOptionalVal(options, "mult", "1.0"));
    summary = bbiTotalSummary(mb->big.bbi);
    high = (int)(summary.maxVal * m);
    low = (int)(summary.minVal * m);
    size = high - low + 1;
    if (size < 1)
	errAbort("need to specify a better low/high");
    AllocArray(counts, size);
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbwList = perBaseWigLoadContinue(mb->big.bbi, section->chrom, section->chromStart, 
							      section->chromEnd);
	struct perBaseWig *pbw;
	for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
	{
	    int len = pbw->chromEnd - pbw->chromStart;
	    for (i = 0; i < len; i++)
		counts[(int)(pbw->data[i]*m) - low]++;
	}
	perBaseWigFreeList(&pbwList);
    }
    output = mustOpen(outputfile, "w");
    j = low;
    for (i = 0; i < size; i++)
	fprintf(output, "%d\t%ld\n", j++, counts[i]);
    carefulClose(&output);
    freeMem(counts);
    metaBigClose(&mb);
}
