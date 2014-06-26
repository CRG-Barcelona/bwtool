/* bwtool_remove - code removing specified parts of a bigWig */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <jkweb/common.h>
#include <jkweb/linefile.h>
#include <jkweb/hash.h>
#include <jkweb/rangeTree.h>
#include <jkweb/options.h>
#include <jkweb/sqlNum.h>
#include <jkweb/basicBed.h>
#include <jkweb/bigWig.h>
#include <beato/bigs.h>
#include "bwtool.h"
#include "bwtool_shared.h"

#define NANUM sqrt(-1)

void usage_remove()
/* Explain usage and exit. */
{
errAbort(
  "bwtool - Data operations on bigWig files\n"
  "usage:\n"
  "   bwtool remove <operator> <value|mask.bed> input.bw[:chr:start-end] output.bw\n"
  "where:\n"
  "   operator is one of the following: \"less\", \"more\", or \"mask\".  In the case\n"
  "   of \"mask\", the next parameter is a bed file containing the regions to remove.\n"
  "   If the operator is not \"mask\", it is a typical binary operator for a thresholding\n"
  "   operation using the value parameter.\n"
  "options:\n"
  "   -inverse   remove the data NOT specified in the operation\n"
  );
}

static struct hash *load_range_tree( char *range_file)
/* load the masking ranges into the tree*/
{
    struct hash *rt_hash = hashNew(0);
    struct bed *beds = bedLoadNAll(range_file, 3);
    struct bed *bed;
    for (bed = beds; bed != NULL; bed = bed->next)
    {
	boolean new_chrom = FALSE;
	struct rbTree *chrom_range = hashFindVal(rt_hash, bed->chrom);
	if (!chrom_range)
	{
	    new_chrom = TRUE;
	    chrom_range = rangeTreeNew();
	}
	rangeTreeAdd(chrom_range, bed->chromStart, bed->chromEnd);
	if (new_chrom)
	    hashAdd(rt_hash, bed->chrom, chrom_range);
    }
    bedFreeList(&beds);
    return rt_hash;
}

static void bwtool_remove_thresh(struct metaBig *mb, enum bw_op_type op, char *val_s,
				 char *outputfile, enum wigOutType wot, unsigned decimals,
                                 boolean condense, boolean wig_only)
/* deal with the thresholding type of removal. */
{
    char wigfile[512];
    safef(wigfile, sizeof(wigfile), "%s.tmp.wig", outputfile);
    FILE *out = mustOpen(wigfile, "w");
    double val = (double)((float)sqlDouble(val_s));
    struct bed *section;
    const double na = NANUM;
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbwList = perBaseWigLoadContinue(mb, section->chrom, section->chromStart, section->chromEnd);
	struct perBaseWig *pbw;
	for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
	{
	    int size = pbw->len;
	    int i;
	    for (i = 0; i < size; i++)
	    {
		switch (op)
		{
		case less:
		{
		    if (pbw->data[i] < val)
			pbw->data[i] = na;
		    break;
		}
		case less_equal:
		{
		    if (pbw->data[i] <= val)
			pbw->data[i] = na;
		    break;
		}
		case more:
		{
		    if (pbw->data[i] > val)
			pbw->data[i] = na;
		    break;
		}
		case more_equal:
		{
		    if (pbw->data[i] >= val)
			pbw->data[i] = na;
		    break;
		}
		case equal:
		{
		    if (pbw->data[i] == val)
			pbw->data[i] = na;
		    break;
		}
		case not_equal:
		{
		    if (pbw->data[i] != val)
			pbw->data[i] = na;
		    break;
		}
		case invalid:
		case mask:
		default:
		{
		    break;
		}
		}
	    }
	}
	perBaseWigOutputNASkip(pbwList, out, wot, decimals, NULL, FALSE, condense);
	perBaseWigFreeList(&pbwList);
    }
    carefulClose(&out);
    if (wig_only)
	rename(wigfile, outputfile);
    else
    {
	writeBw(wigfile, outputfile, mb->chromSizeHash);
	remove(wigfile);
    }
}

static void bwtool_remove_mask(struct metaBig *mb, char *mask_file, char *outputfile, enum wigOutType wot,
			       unsigned decimals, boolean condense, boolean wig_only, boolean inverse)
/* masking */
{
    char wigfile[512];
    safef(wigfile, sizeof(wigfile), "%s.tmp.wig", outputfile);
    FILE *out = mustOpen(wigfile, "w");
    struct hash *rt_hash = load_range_tree(mask_file);
    struct bed *section;
    const double na = NANUM;
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbwList = perBaseWigLoadContinue(mb, section->chrom, section->chromStart, section->chromEnd);
	struct perBaseWig *pbw;
	struct rbTree *chrom_tree = (struct rbTree *)hashFindVal(rt_hash, section->chrom);
	if (chrom_tree && pbwList)
	{
	    for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
	    {
		int size = pbw->chromEnd - pbw->chromStart;
		int i;
		for (i = 0; i < size; i++)
		{
		    /* algorithmically, probably not the greatest solution to do a tree lookup at ever base */
		    /* maybe get back to this another day */
		    boolean does_overlap = rangeTreeOverlaps(chrom_tree, pbw->chromStart + i, pbw->chromStart + i + 1);
		    if ((does_overlap && !inverse) || (!does_overlap && inverse))
			pbw->data[i] = na;
		}
	    }
	    perBaseWigOutputNASkip(pbwList, out, wot, decimals, NULL, FALSE, condense);
	    perBaseWigFreeList(&pbwList);
	}
    }
    carefulClose(&out);
    if (wig_only)
	rename(wigfile, outputfile);
    else
    {
	writeBw(wigfile, outputfile, mb->chromSizeHash);
	remove(wigfile);
    }
    hashFree(&rt_hash);
}

void bwtool_remove(struct hash *options, char *favorites, char *regions, unsigned decimals, enum wigOutType wot,
		   boolean condense, boolean wig_only, char *thresh_type, char *val_or_file, char *bigfile, char *tmp_dir,
		   char *outputfile)
/* bwtool_remove - main for removal program */
{
    boolean inverse = (hashFindVal(options, "inverse") != NULL) ? TRUE : FALSE;
    enum bw_op_type op= get_bw_op_type(thresh_type, inverse);
    if (op == invalid)
	usage_remove();
    struct metaBig *mb = metaBigOpen_check(bigfile, tmp_dir, regions);
    if (op == mask)
	bwtool_remove_mask(mb, val_or_file, outputfile, wot, decimals, condense, wig_only, inverse);
    else
	bwtool_remove_thresh(mb, op, val_or_file, outputfile, wot, decimals, condense, wig_only);
    metaBigClose(&mb);
}

