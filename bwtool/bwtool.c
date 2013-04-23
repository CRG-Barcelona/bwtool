/* bwtool - a bigWig manipulation tool */

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "basicBed.h"
#include "bigWig.h"
#include "bigs.h"
#include "bwtool.h"

void usage()
/* Explain usage and exit. */
{
errAbort(
  "bwtool - Data operations on bigWig files\n"
  "usage:\n"
  "   bwtool command [additional command parameters]\n" 
  "commands:\n"
  "   remove         remove data equal to or thresholded on a given value\n"
  "                  or remove data using ranges specified in a bed file\n"
  "   fill           fill in regions of genome where no data exists with a value\n"
  "   shift          move data on the chromosome\n"
  "   paste          output data from multiple bigWigs and align them one per column\n"
  "                  in tab-delimited output meant to feed into computations\n"
  "   find           find regions of bigWig with given properties\n"
  "   matrix         extract same-sized sections from bigWig to examine as a matrix\n"
  "   summary        provide some summary stats for each region in a bed file\n"
  "   random         print out random data or random regions from the bigWig file\n"
  "   lift           project data from one genome assembly to another using a\n"
  "                  liftOver file (can be lossy)\n"
  "   aggregate      produce plot data as an average of values around the regions\n"
  "                  specified in a bed file\n"
  "   chromgraph     roughly convert to the chromgraph format, suitable for UCSC's\n"
  "                  Genome Graphs page\n"
  "   distribution   produce plot data as the frequency of values seen in the bigWig\n"
  "   sax            run symbolic aggregate approximation (SAX) algorithm on data\n"
  "general options:\n"
  " -wigtype=<bg|fix|var>    output bedGraph, fixedStep, or variableStep wig\n"
  " -regions=bed             use specific regions\n"
  " -condense                condense output, particularly bedGraphs\n"
  " -decimals=n              output specified number of decimals (default 2)\n"
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
struct hash *options = optionParseIntoHash(&argc, argv, FALSE);
/* common options */
char *regions = (char *)hashOptionalVal(options, "regions", NULL);
char *favorites = (char *)hashOptionalVal(options, "favorites", NULL);
unsigned decimals = sqlUnsigned((char *)hashOptionalVal(options, "decimals", "2"));
enum wigOutType wot = get_wig_out_type((char *)hashOptionalVal(options, "wigtype", "fix"));
boolean condense = (hashFindVal(options, "condense") != NULL) ? TRUE : FALSE;
if (argc <= 1)
    usage();
if (sameString(argv[1], "remove"))
{
    if (argc != 6)
	usage_remove();
    else
	bwtool_remove(options, favorites, regions, decimals, wot, condense, argv[2], argv[3], argv[4], argv[5]);
}
else if (sameString(argv[1], "fill"))
{
    if (argc != 5) 
	usage_fill();
    else
	bwtool_fill(options, favorites, regions, decimals, wot, condense, argv[2], argv[3], argv[4]);
}
else if (sameString(argv[1], "shift"))
{
    if (argc != 6) 
	usage_shift();
    else
	bwtool_shift(options, favorites, regions, decimals, wot, condense, argv[3], argv[2], argv[4], argv[5]);
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
	    bwtool_find_extrema(options, favorites, regions, decimals, argv[3], argv[4]);
    }
    else if (sameString(argv[2], "less") || sameString(argv[2], "less-equal") || sameString(argv[2], "more") ||
	     sameString(argv[2], "more-equal") || sameString(argv[2], "equal") || sameString(argv[2], "not"))
    {
	if (argc != 6)
	    usage_find();
	else
	    bwtool_find_thresh(options, favorites, regions, argv[2], argv[3], argv[4], argv[5]);
    }
}
else if (sameString(argv[1], "matrix"))
{
    if (argc != 5)
	usage_matrix();
    else
	bwtool_matrix(options, favorites, argv[2], decimals, argv[3], argv[4]);
}
else if (sameString(argv[1], "distribution") || sameString(argv[1], "dist"))
{
    if (argc != 4)
	usage_distrib();
    else
	bwtool_distrib(options, favorites, regions, decimals, argv[2], argv[3]);
}
/* else if (sameString(argv[1], "autocorr")) */
/* { */
/*     if (argc != 4) */
/* 	usage_autocorr(); */
/*     else */
/* 	bwtool_autocorr(options, favorites, regions, decimals, argv[3], argv[4]); */
/* } */
else if (sameString(argv[1], "random"))
{
    if (argc != 6)
	usage_random();
    else
	bwtool_random(options, favorites, regions, decimals, argv[2], argv[3], argv[4], argv[5]);
}
else if (sameString(argv[1], "aggregate") || sameString(argv[1], "agg"))
{
    if (argc != 6)
	usage_aggregate();
    else
	bwtool_aggregate(options, regions, decimals, argv[2], argv[3], argv[4], argv[5]);
}
else if (sameString(argv[1], "chromgraph") || sameString(argv[1], "cg"))
{
    if (argc != 4)
	usage_chromgraph();
    else
	bwtool_chromgraph(options, favorites, regions, decimals, argv[2], argv[3]);
}
else if (sameString(argv[1], "paste"))
{
    if (argc < 4)
	usage_paste();
    else
    {
	struct slName *list = NULL;
	int i;
	for (i = 3; i < argc; i++)
	{
	    struct slName *name = slNameNew(argv[i]);
	    slAddHead(&list, name);
	}
	slReverse(&list);
	bwtool_paste(options, favorites, argv[2], decimals, list);
    }
}
/* else if (sameString(argv[1], "distance")) */
/* { */
/*     if (argc != 5) */
/* 	usage_distance(); */
/*     else */
/* 	bwtool_distance(options, favorites, regions, decimals, argv[2], argv[3], argv[4]); */
/* } */
else if (sameString(argv[1], "lift"))
{
    if (argc != 5)
	usage_lift();
    else
	bwtool_lift(options, favorites, regions, decimals, wot, argv[2], argv[3], argv[4]);
}
else if (sameString(argv[1], "summary"))
{
    if (argc != 5)
	usage_summary();
    else
	bwtool_summary(options, favorites, regions, decimals, argv[2], argv[3], argv[4]);
}
else if (sameString(argv[1], "sax"))
{
    if (argc != 5)
	usage_sax();
    else
	bwtool_sax(options, favorites, regions, decimals, argv[2], argv[3], argv[4]);
}
else 
    usage();
hashFree(&options);
return 0;
}
