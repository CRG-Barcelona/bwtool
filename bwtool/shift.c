/* bwtool_shift - shifting data  */

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

#include <math.h>

#define NANUM sqrt(-1)

void usage_shift()
/* Explain usage and exit. */
{
errAbort(
  "bwtool shift - move the data on the chromosome by N number of bases\n"
  "   where N can be negative.\n"
  "usage:\n"
  "   bwtool shift <up|down> N input.bw[:chr:start-end] output.bw\n" 
  );
}

void bwtool_shift(struct hash *options, char *favorites, char *regions, unsigned decimals, enum wigOutType wot,
		  boolean condense, char *val_s, char *up_s, char *bigfile, char *outputfile)
/* bwtool_shift - main for shifting program */
{
    const double na = NANUM;
    int shft = sqlSigned(val_s);
    int abs_shft = abs(shft);
    struct metaBig *mb = metaBigOpen_favs(bigfile, regions, favorites);
    if (!mb)
	errAbort("problem opening %s", bigfile);
    char wigfile[512];
    safef(wigfile, sizeof(wigfile), "%s.tmp.wig", outputfile);
    FILE *out = mustOpen(wigfile, "w");
    struct bed *section;
    boolean up = TRUE;
    if (up_s && sameString("up", up_s))
	up = TRUE;
    else if (up_s && sameString("down", up_s))
	up = FALSE;
    else
	errAbort("must choose \"up\" or \"down\" with shift function");
    if (shft == 0)
	errAbort("it doesn't make sense to shift by zero.");
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, section->chrom, section->chromStart, 
							      section->chromEnd, FALSE, na);
	int i;
	/* if the shift size is bigger than the section, NA the entire thing */
	int size = pbw->len;
	if (abs_shft >= size)
	    for (i = 0; i < size; i++)
		pbw->data[i] = na;
	else
	{
	    if (!up)
	    {
		for (i = size-1; i >= abs_shft; i--)
		    pbw->data[i] = pbw->data[i - abs_shft];
		for (; i >= 0; i--)
		    pbw->data[i] = na;
	    }
	    else
	    {
		for (i = 0; i < size - abs_shft; i++)
		    pbw->data[i] = pbw->data[i + abs_shft];
		for (; i < size; i++)
		    pbw->data[i] = na;
	    }
	}
	perBaseWigOutputNASkip(pbw, out, wot, decimals, NULL, FALSE, condense);
	perBaseWigFree(&pbw);
    }
    carefulClose(&out);
    writeBw(wigfile, outputfile, mb->chromSizeHash);
    remove(wigfile);
    metaBigClose(&mb);
}
