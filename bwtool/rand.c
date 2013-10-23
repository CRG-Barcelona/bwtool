/* bwtool_random - retrieve random data */

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "basicBed.h"
#include "bigWig.h"
#include "bigs.h"
#include "random_coord.h"
#include "bwtool.h"

#include <stdlib.h>

void usage_random()
/* Explain usage of randomized data program and exit. */
{
errAbort(
  "bwtool random - retrieve data from random regions in the bigWig.\n"
  "usage:\n"
  "   bwtool random n size input.bw output.txt\n" 
  "where n is the number of regions desired, and size is the size of each\n"
  "region retrieved.\n"
  "options:\n"
  "   -bed            instead of outputting data, just output the bed indicating\n"
  "                   where the data lies.\n"
  "   -NA-perc=p      maximum percent of region allowed to contain NA values\n"
  "   -blacklist=bed  specifically avoid these regions."
  );
}

void bwtool_random(struct hash *options, char *favorites, char *regions, unsigned decimals, 
		   double fill, char *num_s, char *size_s, char *bigfile, char *output_file)
/* random - main ... random number generation takes place here.  */
{
    struct metaBig *mb = metaBigOpen_favs(bigfile, regions, favorites);
    FILE *out = mustOpen(output_file, "w");
    boolean just_bed = (hashFindVal(options, "bed") != NULL) ? TRUE : FALSE;
    unsigned N = sqlUnsigned(num_s);
    unsigned size = sqlUnsigned(size_s);
    char *blacklist_file = hashFindVal(options, "blacklist");
    struct bed *blacklist = NULL;
    struct perBaseWig *pbw;
    struct perBaseWig *pbwList = NULL;
    double NA_perc = sqlDouble(hashOptionalVal(options, "NA-perc", "0.4"));
    struct random_coord *rc = NULL;
    unsigned long max = 0;
    unsigned kept = 0;
    srand(time(NULL));
    if (blacklist_file)
	blacklist = bedLoadNAll(blacklist_file, 3);
    rc = random_coord_init(mb->chromSizeHash, blacklist);
    max = rc->length - size;
    while (kept < N)
    {
	struct bed *bed = NULL;
	while (bed == NULL)
	{
	    unsigned long rand = random_in_range(0, max);
	    bed = random_bed(rc, size, rand);
	}
	pbw = perBaseWigLoadSingleContinue(mb, bed->chrom, bed->chromStart, bed->chromEnd, FALSE, fill);
	int size_pbw = pbw->chromEnd - pbw->chromStart;
	int count_NA = 0;
	int i;
	if (isnan(fill))
	{
	    for (i = 0; i < size_pbw; i++)
	    {
		if (isnan(pbw->data[i]))
		    count_NA++;
	    }
	    if ((double)count_NA/size_pbw <= NA_perc)
	    {
		slAddHead(&pbwList, pbw);
		kept++;
	    }
	}
	else
	{
	    slAddHead(&pbwList, pbw);
	    kept++;
	}	    
    }
    slSort(&pbwList, bedCmp);
    /* output either the bed or a tab-separated list of vals*/ 
    for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
    {
	if (just_bed)
	    fprintf(out, "%s\t%d\t%d\n", pbw->chrom, pbw->chromStart, pbw->chromEnd);
	else
	{
	    int i;
	    for (i = 0; i < (pbw->chromEnd - pbw->chromStart); i++)
		fprintf(out, "%0.*f%c", decimals, pbw->data[i], (i < pbw->chromEnd - pbw->chromStart - 1) ? '\t' : '\n');
	}
    }
    carefulClose(&out);
    bedFreeList(&blacklist);
    random_coord_free(&rc);
}
