/* bwtool_split - splitting the bigWig into useful pieces */

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

void usage_split()
/* Explain usage of the splitting program and exit. */
{
errAbort(
  "bwtool split - split the bigWig into chunks\n"
  "usage:\n"
  "   bwtool split size file.bw[:chr:start-end] prefix\n"
  "options:\n"
  "   -min-gap=g     if gap is smaller than this, continue to add it\n"
  "                  to the current piece (default 1)\n"
  );
}

struct bed *newBed(char *chrom, int start, int end)
{
    struct bed *bed;
    AllocVar(bed);
    bed->chrom = cloneString(chrom);
    bed->chromStart = start;
    bed->chromEnd = end;
    return bed;
}

void bwtool_split(struct hash *options, char *regions, char *size_s, char *bigfile, char *tmp_dir, char *outputfile)
/* bwtool_split - main for the splitting program */
{
    struct metaBig *mb = metaBigOpenWithTmpDir(bigfile, tmp_dir, regions);
    FILE *output = mustOpen(outputfile, "w");
    struct bed *section;
    struct bed *splitList = NULL;
    int size = 0;
    unsigned min_gap = sqlUnsigned((char *)hashOptionalVal(options, "min_gap", "1"));
    unsigned chunk_size = sqlUnsigned(size_s);
    char chrom[256] = "";
    int start = -1, end = 0;
    boolean over_size = FALSE;
    int ix = 1;
    int gap = 0;
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbwList = perBaseWigLoadContinue(mb, section->chrom, section->chromStart,
							      section->chromEnd);
	struct perBaseWig *pbw;
	for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
	{
	    int length = pbw->chromEnd - pbw->chromStart;
	    if (end > 0)
		gap = pbw->chromStart - end;
	    if (!sameString(chrom, pbw->chrom))
	    {
		if (!sameString(chrom, ""))
		    slAddHead(&splitList, newBed(chrom, start, end));
		strcpy(chrom, pbw->chrom);
		start = pbw->chromStart;
		end = pbw->chromEnd;
		if (size + length > chunk_size)
		    size = length;
		else
		    size += length;
	    }
	    else
	    {
		if ((size + length + gap > chunk_size) && (gap >= min_gap))
		{
		    slAddHead(&splitList, newBed(chrom, start, end));
		    start = pbw->chromStart;
		    end = pbw->chromEnd;
		    size = length;
		}
		else
		{
		    size += length + gap;
		    end = pbw->chromEnd;
		}
	    }
	}
	perBaseWigFreeList(&pbwList);
    }
    slAddHead(&splitList, newBed(chrom, start, end));
    slReverse(&splitList);
    for (section = splitList; section != NULL; section = section->next)
    {
	fprintf(output, "%s\t%d\t%d\n", section->chrom, section->chromStart, section->chromEnd);
    }
    carefulClose(&output);
    metaBigClose(&mb);
    bedFreeList(&splitList);
}
