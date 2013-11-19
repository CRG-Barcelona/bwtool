/* bwtool_fill - code adding filler data to  */

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

void usage_fill()
/* Explain usage and exit. */
{
errAbort(
  "bwtool fill - fill given sections or whole chromosomes with\n"
  "   a given value anywhere there is no data.\n"
  "usage:\n"
  "   bwtool fill <val> input.bw[:chr:start-end] output.bw\n" 
  );
}

void bwtool_fill(struct hash *options, char *favorites, char *regions, unsigned decimals, enum wigOutType wot,
		 boolean condense, char *val_s, char *bigfile, char *outputfile)
/* bwtool_fill - main for filling program */
{
    double val = sqlDouble(val_s);
    struct metaBig *mb = metaBigOpen_check(bigfile, regions);
    char wigfile[512];
    safef(wigfile, sizeof(wigfile), "%s.tmp.wig", outputfile);
    FILE *out = mustOpen(wigfile, "w");
    struct bed *section;
    int i;
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, section->chrom, section->chromStart, 
							      section->chromEnd, FALSE, val);
	perBaseWigOutput(pbw, out, wot, decimals, NULL, FALSE, condense);
	perBaseWigFree(&pbw);
    }
    carefulClose(&out);
    writeBw(wigfile, outputfile, mb->chromSizeHash);
    remove(wigfile);
    metaBigClose(&mb);
}
