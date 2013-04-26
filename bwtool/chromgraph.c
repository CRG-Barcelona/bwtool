/* bwtool_chromgraph - chromgraph plotting  */

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "basicBed.h"
#include "bigWig.h"
#include "bigs.h"
#include "bwtool.h"

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
		    char *bigfile, char *outputfile)
/* bwtool_chromgraph - main for making the chromgraph file */
{
    struct metaBig *mb = metaBigOpen_favs(bigfile, regions, favorites);
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
	struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, section->chrom, section->chromStart, section->chromEnd, FALSE);

	windowPos = pbw->chromStart;
	while (windowPos < pbw->chromEnd)
	{
	    numBases = 0;
	    sum = 0;
	    for (i = windowPos; (i < pbw->chromEnd) && (i < windowPos + every); i++)
	    {
		if (!isnan(pbw->data[i]))
		{
		    numBases++;
		    sum += pbw->data[i];
		}
	    }
	    if (numBases > 0)
		sum = sum / numBases;
	    fprintf(output, "%s\t%d\t%0.*f\n", pbw->chrom, windowPos, decimals, sum);
	    windowPos += every;
	} 
	perBaseWigFreeList(&pbw);
    }
    carefulClose(&output);
    metaBigClose(&mb);
}
