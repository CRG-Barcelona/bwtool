/* wigmake - bigWig -> wig... shouldn't be run directly */
/*   ** this isn't a full-featured program.  I'ts meant */
/*   ** just for the test script. */

/* It's run like */
/*   wigmake bw wig decimals wigtype condense */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif 

#include "common.h"
#include "basicBed.h"
#include "hash.h"
#include "options.h"
#include "linefile.h"
#include "sqlNum.h"
#include "bigWig.h"
#include "bwgInternal.h"
#include "bigs.h"
#include "metaBig.h"

#define NANUM sqrt(-1)

int main(int argc, char *argv[])
/* Process command line. */
{
    if (argc != 6)
	errAbort("bad running of bwmake");
    struct metaBig *mb = metaBigOpen(argv[1], NULL);
    FILE *out = mustOpen(argv[2], "w");
    int decimals = sqlSigned(argv[3]);
    struct bed *section;
    double na = NANUM;
    boolean condense = FALSE;
    enum wigOutType wot = get_wig_out_type(argv[4]);
    if (sameString(argv[5], "condense"))
	condense = TRUE;
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, section->chrom, section->chromStart, 
							      section->chromEnd, FALSE, na);
	perBaseWigOutputNASkip(pbw, out, wot, decimals, NULL, FALSE, condense);
	perBaseWigFree(&pbw);
    }
    carefulClose(&out);
    metaBigClose(&mb);
    return 0;
}
