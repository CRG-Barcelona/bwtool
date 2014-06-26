/* bwtool_roll - a special rolling-mean,etc case of window */

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

void usage_roll()
/* Explain usage of the rolling-average program and exit */
{
errAbort(
  "bwtool roll - slide a window across the bigWig and print a rolling mean\n"
  "   or total:\n"
  "usage:\n"
  "   bwtool roll <command> size file.bw output.txt\n"
  "where command is either \"mean\" or \"total\""
  "options:\n"
  "   -max-NA       maximum NA-valued bases to consider a region legitimate.\n"
  "   -min-mean=m   remove regions in output having calculated means < m\n"
  );
}

void add_to_tots(double d, int *p_num_na, double *p_total)
/* add to front  */
{
    if (isnan(d))
	*p_num_na = *p_num_na+ 1;
    else
	*p_total = *p_total + d;
}

void sub_from_tots(double d, int *p_num_na, double *p_total)
/* take off the back */
{
    if (isnan(d))
	*p_num_na = *p_num_na - 1;
    else *p_total = *p_total - d;
}

enum roll_command
{
    roll_mean = 0,
    roll_total = 1,
};

void bwtool_roll(struct hash *options, char *favorites, char *regions, unsigned decimals, double fill,
		 enum wigOutType wot, char *command, char *size_s, char *bigfile, char *tmp_dir, char *outputfile)
/* bwtool_roll - main for the rolling-mean program */
/* this function is too long. it'd be nice to break it up some time. */
{
    struct metaBig *mb = metaBigOpenWithTmpDir(bigfile, tmp_dir, regions);
    int step = (int)sqlUnsigned((char *)hashOptionalVal(options, "step", "1"));
    int max_na = (int)sqlSigned((char *)hashOptionalVal(options, "max-NA", "-1"));
    char *min_mean_s = (char *)hashOptionalVal(options, "min-mean", "unused");
    double min_mean = -DBL_MAX;
    if (!sameString(min_mean_s,"unused"))
	min_mean = sqlDouble(min_mean_s);
    int size = sqlSigned(size_s);
    if (max_na == -1)
	max_na = size/2;
    else if (max_na > size)
	max_na = size - 1;
    if (size < 1)
	errAbort("size must be >= 1 for bwtool window");
    FILE *out = (outputfile) ? mustOpen(outputfile, "w") : stdout;
    struct bed *section;
    boolean broken = TRUE;  /* for headers */
    enum roll_command com;
    if (sameWord(command, "mean"))
	com = roll_mean;
    else if (sameWord(command, "total"))
	com = roll_total;
    else
	errAbort("Pick a roll command: mean or total");
    for (section = mb->sections; section != NULL; section = section->next)
    {
	if (size <= section->chromEnd - section->chromStart)
	{
	    struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, section->chrom, section->chromStart,
								  section->chromEnd, FALSE, fill);
	    int i, j;
	    int len = section->chromEnd - section->chromStart;
	    double total = 0;
	    int num_na = 0;
	    /* load data */
	    for (i = 0; i < size; i++)
		add_to_tots(pbw->data[i], &num_na, &total);
	    i = 0;
	    do
	    {
		int s = pbw->chromStart + i;
		int e = pbw->chromStart + i + size;
 		int st = step;
		double mean = total/(size - num_na);
		/* the next two calculations center it */
		s += size/2 - step/2;
		e = s + step;
		/* output */
		if ((num_na <= max_na) && (mean >= min_mean))
		{
		    double out_val;
		    if (com == roll_mean)
			out_val = mean;
		    else
			out_val = total;
		    if (wot == fixStepOut)
		    {
			if (broken)
			    fprintf(out, "fixedStep chrom=%s start=%d step=%d span=%d\n", pbw->chrom, s+1, step, step);
			fprintf(out, "%0.*f\n", decimals, out_val);
		    }
		    else if (wot == varStepOut)
		    {
			if (broken)
			    fprintf(out, "variableStep chrom=%s span=%d\n", pbw->chrom, step);
			fprintf(out, "%d\t%0.*f\n", s+1, decimals, out_val);
		    }
		    else
			fprintf(out, "%s\t%d\t%d\t%0.*f\n", pbw->chrom, s, e, decimals, out_val);
		    broken = FALSE;
		}
		else
		    broken = TRUE;
		/* move */
		while (st > 0)
		{
		    if (i + size <= pbw->len)
		    {
			add_to_tots(pbw->data[i+size], &num_na, &total);
			sub_from_tots(pbw->data[i], &num_na, &total);
		    }
		    else
			break;
		    i++;
		    st--;
		}
	    } while (i <= pbw->len - size);
	    perBaseWigFree(&pbw);
	}
    }
    metaBigClose(&mb);
    carefulClose(&out);
}
