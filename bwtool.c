/* bwtool - a bigWig manipulation tool */

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

#define NANUM sqrt(-1)

void usage()
/* Explain usage and exit. */
{
errAbort(
  "bwtool - Data operations on bigWig files\n"
  "usage:\n"
  "   bwtool command [additional command parameters]\n"
  "commands:\n"
  "   aggregate      (or \"agg\") produce plot data as an average of values around\n"
  "                  the regions specified in a bed file\n"
  "   chromgraph     roughly convert to the chromgraph format, suitable for UCSC's\n"
  "                  Genome Graphs page\n"
  "   distribution   (or \"dist\") produce plot data as the frequency of values seen\n"
  "                  in the bigWig\n"
  "   extract        (or \"ex\") extract data in some other ways than paste, matrix, or\n"
  "                  window with a given bed, preserving strand directionality.\n"
  "   fill           fill in regions of genome where no data exists with a value\n"
  "   find           find regions of bigWig with given properties\n"
  "   lift           project data from one genome assembly to another using a\n"
  "                  liftOver file (can be lossy)\n"
  "   matrix         extract same-sized sections from bigWig to examine as a matrix\n"
  "   paste          output data from multiple bigWigs and align them one per column\n"
  "                  in tab-delimited output meant to feed into computations\n"
/* #ifdef HAVE_LIBGSL */
/*   "   random         print out random data or random regions from the bigWig file\n" */
/* #endif */
  "   remove         remove data equal to or thresholded on a given value\n"
  "                  or remove data using ranges specified in a bed file\n"
  "   roll           compute rolling means, etc\n"
  "   sax            run symbolic aggregate approximation (SAX) algorithm on data\n"
  "   shift          move data on the chromosome\n"
/*   "   split          make a set of files describing evenly-sized regions of the bigWig,\n" */
/*   "                  each of which may be used separately on a cluster and combined later\n" */
  "   summary        provide some summary stats for each region in a bed file\n"
  "   window         print out tiling windows of data in comma-separated lists\n\n"
  "general options:\n"
  " -wigtype=<bg|fix|var>    output bedGraph, fixedStep, or variableStep wig\n"
  " -wig-only                for bigWig-creating programs, make a wig instead\n"
  " -regions=bed             use specific regions\n"
  " -condense                condense output, particularly bedGraphs\n"
  " -decimals=n              output specified number of decimals (default 2)\n"
  " -fill=val                some programs allow filling missing parts of the bigWig\n"
  "                          with a specified value prior to using data.\n"
  " -pseudo=val              add a pseudo-count at every value\n"
  " -o=output.txt            where normally standard output is written, write to a\n"
  "                          file instead.\n"
  " -tmp-dir=dir             by default, bigWig caching is done in /tmp/udcCache/*.\n"
  "                          Override this by setting dir to the desired path.\n"
  );
}

enum bw_op_type get_bw_op_type(char *thresh_type, boolean inverse)
/* speed things up inside loop.  we don't want a string comparison for every datapoint */
/* to know what kind of operation to perform. */
{
    if (thresh_type)
    {
	if (sameWord(thresh_type, "less"))
	{
	    if (inverse)
		return more_equal;
	    return less;
	}
	else if (sameWord(thresh_type, "less-equal"))
	{
	    if (inverse)
		return more;
	    return less_equal;
	}
	else if (sameWord(thresh_type, "more"))
	{
	    if (inverse)
		return less_equal;
	    return more;
	}
	else if (sameWord(thresh_type, "more-equal"))
	{
	    if (inverse)
		return less;
	    return more_equal;
	}
	else if (sameWord(thresh_type, "equal"))
	{
	    if (inverse)
		return not_equal;
	    return equal;
	}
	else if (sameWord(thresh_type, "not-equal"))
	{
	    if (inverse)
		return equal;
	    return not_equal;
	}
	else if (sameWord(thresh_type, "mask"))
	    return mask;
    }
    return invalid;
}

int main(int argc, char *argv[])
/* Process command line. */
{
struct hash *options = optionParseIntoHashExceptNumbers(&argc, argv, FALSE);
/* common options */
char *tmp_dir = (char *)hashOptionalVal(options, "tmp-dir", NULL);
char *regions = (char *)hashOptionalVal(options, "regions", NULL);
char *output_file = (char *)hashOptionalVal(options, "o", NULL);
char *favorites = (char *)hashOptionalVal(options, "favorites", NULL);
unsigned decimals = sqlUnsigned((char *)hashOptionalVal(options, "decimals", "2"));
enum wigOutType wot = get_wig_out_type((char *)hashOptionalVal(options, "wigtype", "fix"));
boolean condense = (hashFindVal(options, "condense") != NULL) ? TRUE : FALSE;
boolean wig_only = (hashFindVal(options, "wig-only") != NULL) ? TRUE : FALSE;
char *fill_s = hashFindVal(options, "fill");
double na = NANUM;
double fill = na;
if (fill_s)
    fill = sqlDouble(hashFindVal(options, "fill"));
if (argc <= 1)
    usage();
if (sameString(argv[1], "remove"))
{
    if (argc != 6)
	usage_remove();
    else
	bwtool_remove(options, favorites, regions, decimals, wot, condense, wig_only, argv[2], argv[3], argv[4], tmp_dir, argv[5]);
}
else if (sameString(argv[1], "fill"))
{
    if (argc != 5)
	usage_fill();
    else
	bwtool_fill(options, favorites, regions, decimals, wot, condense, argv[2], argv[3], tmp_dir, argv[4]);
}
else if (sameString(argv[1], "shift"))
{
    if (argc != 5)
	usage_shift();
    else
	bwtool_shift(options, favorites, regions, decimals, wot, condense, argv[2], argv[3], tmp_dir, argv[4]);
}
else if (sameString(argv[1], "find"))
{
    if (argc < 5)
	usage_find();
    else if (sameString(argv[2], "local-extrema"))
    {
	if (argc != 5)
	    usage_find();
	else
	    bwtool_find_extrema(options, favorites, regions, decimals, fill, argv[3], tmp_dir, argv[4]);
    }
    else if (sameString(argv[2], "less") || sameString(argv[2], "less-equal") || sameString(argv[2], "more") ||
	     sameString(argv[2], "more-equal") || sameString(argv[2], "equal") || sameString(argv[2], "not"))
    {
	if (argc != 6)
	    usage_find();
	else
	    bwtool_find_thresh(options, favorites, regions, fill, argv[2], argv[3], argv[4], tmp_dir, argv[5]);
    }
    else if (sameString(argv[2], "maxima"))
    {
	if (argc != 6)
	    usage_find();
	else
	    bwtool_find_max(options, favorites, argv[3], fill, argv[4], tmp_dir, argv[5]);
    }
}
else if (sameString(argv[1], "matrix"))
{
    if (argc != 6)
	usage_matrix();
    else
	bwtool_matrix(options, favorites, argv[3], decimals, fill, argv[2], argv[4], tmp_dir, argv[5]);
}
else if (sameString(argv[1], "distribution") || sameString(argv[1], "dist"))
{
    if (argc != 4)
	usage_distrib();
    else
	bwtool_distrib(options, favorites, regions, decimals, argv[2], tmp_dir, argv[3]);
}
/* #ifdef HAVE_LIBGSL  */
/* else if (sameString(argv[1], "random")) */
/* { */
/*     if (argc != 6) */
/* 	usage_random(); */
/*     else */
/* 	bwtool_random(options, favorites, regions, decimals, fill, argv[2], argv[3], argv[4], argv[5]); */
/* } */
/* #endif */
else if (sameString(argv[1], "aggregate") || sameString(argv[1], "agg"))
{
    decimals = sqlUnsigned((char *)hashOptionalVal(options, "decimals", "6"));
    if (argc != 6)
	usage_aggregate();
    else
	bwtool_aggregate(options, regions, decimals, fill, argv[2], argv[3], argv[4], tmp_dir, argv[5]);
}
else if (sameString(argv[1], "chromgraph") || sameString(argv[1], "cg"))
{
    if (argc != 4)
	usage_chromgraph();
    else
	bwtool_chromgraph(options, favorites, regions, decimals, fill, argv[2], tmp_dir, argv[3]);
}
else if (sameString(argv[1], "paste"))
{
    if (hashFindVal(options, "wigtype") == NULL)
	wot = bedGraphOut;
    if (argc < 3)
	usage_paste();
    else
    {
	struct slName *list = NULL;
	int i;
	for (i = 2; i < argc; i++)
	{
	    struct slName *name = slNameNew(argv[i]);
	    slAddHead(&list, name);
	}
	slReverse(&list);
	bwtool_paste(options, favorites, regions, decimals, fill, wot, &list, tmp_dir, output_file);
    }
}
else if (sameString(argv[1], "lift"))
{
    if (argc != 5)
	usage_lift();
    else
	bwtool_lift(options, favorites, regions, decimals, wot, argv[2], tmp_dir, argv[3], argv[4]);
}
else if (sameString(argv[1], "roll"))
{
    if (argc != 6)
	usage_roll();
    else
	bwtool_roll(options, favorites, regions, decimals, fill, wot, argv[2], argv[3], argv[4], tmp_dir, argv[5]);
}
else if (sameString(argv[1], "summary"))
{
    if (argc != 5)
	usage_summary();
    else
	bwtool_summary(options, favorites, regions, decimals, fill, argv[2], argv[3], tmp_dir, argv[4]);
}
else if (sameString(argv[1], "sax"))
{
    if (argc != 5)
	usage_sax();
    else
	bwtool_sax(options, favorites, regions, decimals, argv[2], argv[3], tmp_dir, argv[4]);
}
else if (sameString(argv[1], "split"))
{
    if (argc != 5)
	usage_split();
    else
	bwtool_split(options, regions, argv[2], argv[3], tmp_dir, argv[4]);
}
else if (sameString(argv[1], "window") || sameString(argv[1], "win"))
{
    if (argc != 4)
	usage_window();
    else
	bwtool_window(options, favorites, regions, decimals, fill, argv[2], argv[3], tmp_dir, output_file);
}
else if (sameString(argv[1], "extract") || sameString(argv[1], "ex"))
{
    if (argc != 6)
	usage_extract();
    else
	bwtool_extract(options, argv[3], decimals, fill, argv[2], argv[4], tmp_dir, argv[5]);
}
else
    usage();
hashFree(&options);
return 0;
}
