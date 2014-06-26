/* bwtool_paste - simultaneously output same regions of multiple bigWigs */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <jkweb/common.h>
#include <jkweb/obscure.h>
#include <jkweb/linefile.h>
#include <jkweb/hash.h>
#include <jkweb/options.h>
#include <jkweb/sqlNum.h>
#include <jkweb/basicBed.h>
#include <jkweb/bigWig.h>
#include <beato/bigs.h>
#include "bwtool.h"
#include "bwtool_shared.h"
#include <beato/cluster.h>

void usage_window()
/* Explain usage of the window-tiling program and exit */
{
errAbort(
  "bwtool window - slide a window across the bigWig and at each step print data\n"
  "   in a format like:\n"
  "\n"
  "   chrom<TAB>start<TAB>end<TAB>val_start,val_start+1,val_start+2,...,val_end\n"
  "\n"
  "usage:\n"
  "   bwtool window size file.bw\n"
  "options:\n"
  "   -step=n         skip n bases when sliding window (default 1)\n"
  "   -skip-NA        don't output lines (windows) containing any NA values\n"
  "   -center         print start and end coordinates of the middle of the window\n"
  "                   with size step such that the start/ends are connected each\n"
  "                   line (if step < size)\n"
  );
}

void bwtool_window(struct hash *options, char *favorites, char *regions, unsigned decimals,
                   double fill, char *size_s, char *bigfile, char *tmp_dir, char *output_file)
/* bwtool_window - main for the windowing program */
{
    struct metaBig *mb = metaBigOpenWithTmpDir(bigfile, tmp_dir, regions);
    boolean skip_na = (hashFindVal(options, "skip-NA") != NULL) ? TRUE : FALSE;
    if (!isnan(fill) && skip_na)
	errAbort("cannot use -skip_na with -fill");
    boolean center = (hashFindVal(options, "center") != NULL) ? TRUE : FALSE;
    int step = (int)sqlUnsigned((char *)hashOptionalVal(options, "step", "1"));
    int size = sqlSigned(size_s);
    if (size < 1)
	errAbort("size must be >= 1 for bwtool window");
    FILE *out = (output_file) ? mustOpen(output_file, "w") : stdout;
    struct bed *section;
    for (section = mb->sections; section != NULL; section = section->next)
    {
	if (size <= section->chromEnd - section->chromStart)
	{
	    /* I think window should be split into to two functions */
	    /* when skipping NA, perBaseWigLoadContinue should be used */
	    struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, section->chrom, section->chromStart,
								  section->chromEnd, FALSE, fill);
	    int i, j;
	    for (i = 0; i <= pbw->len - size; i += step)
	    {
		int s = pbw->chromStart + i;
		int e = pbw->chromStart + i + size;
		if (center)
		{
		    s += size/2 - step/2;
		    e = s + step;
		}
		boolean has_NA = FALSE;
		if (skip_na)
		{
		    /* stupidly inefficient */
		    for (j = i; j < i + size; j++)
			if (isnan(pbw->data[j]))
			{
			    i = j-step+1;
			    has_NA = TRUE;
			    break;
			}
		}
		if (!has_NA)
 		{
		    fprintf(out, "%s\t%d\t%d\t", pbw->chrom, s, e);
		    for (j = i; j < i + size; j++)
			if (isnan(pbw->data[j]) && (j == i + size - 1))
			    fprintf(out, "NA\n");
			else if (isnan(pbw->data[j]))
			    fprintf(out, "NA,");
			else if (j == i + size - 1)
			    fprintf(out, "%0.*f\n", decimals, pbw->data[j]);
			else
			    fprintf(out, "%0.*f,", decimals, pbw->data[j]);
		}
	    }
	    perBaseWigFree(&pbw);
	}
    }
    metaBigClose(&mb);
    carefulClose(&out);
}
