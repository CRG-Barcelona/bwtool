/* bwtool_chromgraph - chromgraph plotting  */

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

void usage_chromgraph()
/* Explain usage of chromgraph program and exit. */
{
errAbort(
  "bwtool chromgraph - produce a chromgraph file usable by the UCSC Genome Graphs tool\n"
  "usage:\n"
  "   bwtool chromgraph input.bw[:chr:start-end] output.txt\n"
  "options:\n"
  "   -every=N     Output datapoints every N bases instead of default (10,000 bp)\n"
  );
}

void bwtool_chromgraph(struct hash *options, char *favorites, char *regions, unsigned decimals,
		       double fill, char *bigfile, char *tmp_dir, char *outputfile)
/* bwtool_chromgraph - main for making the chromgraph file */
{
    struct metaBig *mb = metaBigOpen_check(bigfile, tmp_dir, regions);
    FILE *output = mustOpen(outputfile, "w");
    struct bed *section;
    struct bbiSummaryElement summary;
    unsigned every = sqlUnsigned((char *)hashOptionalVal(options, "every", "10000"));
    if (mb->type != isaBigWig)
	errAbort("file not bigWig type");
    /* max - min will be the size of the count array */
    summary = bbiTotalSummary(mb->big.bbi);
    for (section = mb->sections; section != NULL; section = section->next)
    {
	int i;
	int numBases = 0;
	double sum = 0;
	int windowPos;
	struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, section->chrom, section->chromStart, section->chromEnd, FALSE, fill);

	windowPos = pbw->chromStart;
	while (windowPos < pbw->chromEnd)
	{
	    numBases = 0;
	    sum = 0;
	    int end = (windowPos + every > pbw->chromEnd) ? pbw->chromEnd : windowPos + every;
	    int middle = windowPos + (end-windowPos)/2;
	    for (i = windowPos; i < end; i++)
	    {
		if (!isnan(pbw->data[i]))
		{
		    numBases++;
		    sum += pbw->data[i];
		}
	    }
	    if (numBases > 0)
		sum = sum / numBases;
	    fprintf(output, "%s\t%d\t%0.*f\n", pbw->chrom, middle, decimals, sum);
	    windowPos += every;
	}
	perBaseWigFreeList(&pbw);
    }
    carefulClose(&output);
    metaBigClose(&mb);
}
