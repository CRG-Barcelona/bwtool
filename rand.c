/* bwtool_random - retrieve random data */
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
#include <beato/random_coord.h>
#include "bwtool.h"
#include "bwtool_shared.h"

#include <stdlib.h>

#ifdef HAVE_LIBGSL

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
  "   -blacklist=bed  specifically avoid these regions.\n"
  "   -seed=s         seed the random number generator with some positive integer\n"
  );
}

void bwtool_random(struct hash *options, char *favorites, char *regions, unsigned decimals,
		   double fill, char *num_s, char *size_s, char *bigfile, char *tmp_dir, char *output_file)
/* random - main ... random number generation takes place here.  */
{
    struct metaBig *mb = metaBigOpen_check(bigfile, tmp_dir, regions);
    FILE *out = mustOpen(output_file, "w");
    boolean just_bed = (hashFindVal(options, "bed") != NULL) ? TRUE : FALSE;
    unsigned seed = sqlUnsigned((char *)hashOptionalVal(options, "seed", "0"));
    unsigned N = sqlUnsigned(num_s);
    unsigned size = sqlUnsigned(size_s);
    char *blacklist_file = hashFindVal(options, "blacklist");
    struct bed *blacklist = NULL;
    struct perBaseWig *pbw;
    struct perBaseWig *pbwList = NULL;
    double NA_perc = sqlDouble(hashOptionalVal(options, "NA-perc", "0.4"));
    if (blacklist_file)
	blacklist = bedLoadNAll(blacklist_file, 3);
    pbwList = random_pbw_list(size, N, mb, NA_perc, fill, blacklist, seed);
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
}

#else

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
  "   -blacklist=bed  specifically avoid these regions.\n"
  "   -seed=s         seed the random number generator with some positive integer\n"
  );
}

void bwtool_random(struct hash *options, char *favorites, char *regions, unsigned decimals,
		   double fill, char *num_s, char *size_s, char *bigfile, char *tmp_dir, char *output_file)
/* random - main ... random number generation takes place here.  */
{
}

#endif
