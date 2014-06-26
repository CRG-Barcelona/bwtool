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
#include <beato/cluster.h>
#include "bwtool_shared.h"

void usage_paste()
/* Explain usage of paste program and exit. */
{
errAbort(
  "bwtool paste - simultaneously output same regions of multiple bigWigs\n"
  "usage:\n"
  "   bwtool paste input1.bw input2.bw input3.bw ...\n"
  "options:\n"
  "   -header           put header with labels from file or filenames\n"
  "   -consts=c1,c2...  add constants to output lines\n"
  "   -consts-means     add means of input bigWigs as constants every line in the\n"
  "                     same order as the bigWigs.\n"
  "   -consts-totals    like -consts-means except the total sum\n"
  "   -consts-covs      like -consts-means except number of bases covered\n"
  "   -skip-NA          don't output lines (bases) where one of the inputs is NA\n"
  "   -skip-min=m       skip in output where one of the inputs is < m\n"
  "                     Care should be taken with this option so that the decimal\n"
  "                     precision corresponds to a threshold.  E.g. if two decimal\n"
  "                     places is the output (the default), then instead of using\n"
  "                     -min=0, use -min=0.01\n"
  );
}

void print_line(struct perBaseWig *pbw_list, struct slDouble *c_list, int decimals, enum wigOutType wot, int i, FILE *out)
{
    struct perBaseWig *pbw;
    struct slDouble *c;
    if (wot == bedGraphOut)
	fprintf(out, "%s\t%d\t%d\t", pbw_list->chrom, pbw_list->chromStart+i, pbw_list->chromStart+i+1);
    else if (wot == varStepOut)
	fprintf(out, "%d\t", pbw_list->chromStart+i+1);
    for (pbw = pbw_list; pbw != NULL; pbw = pbw->next)
    {
	if (isnan(pbw->data[i]))
	    fprintf(out, "NA");
	else
	    fprintf(out, "%0.*f", decimals, pbw->data[i]);
	fprintf(out, "%c", (c_list == NULL) && (pbw->next == NULL) ? '\n' : '\t');
    }
    for (c = c_list; c != NULL; c = c->next)
	fprintf(out, "%0.*f%c", decimals, c->val, (c->next == NULL) ? '\n' : '\t');
}

boolean has_na(struct perBaseWig *pbw_list, int i)
{
    struct perBaseWig *pbw;
    for (pbw = pbw_list; pbw != NULL; pbw = pbw->next)
	if (isnan(pbw->data[i]))
	    return TRUE;
    return FALSE;
}

boolean has_under(struct perBaseWig *pbw_list, int i, double m)
{
    struct perBaseWig *pbw;
    for (pbw = pbw_list; pbw != NULL; pbw = pbw->next)
	if (pbw->data[i] < m)
	    return TRUE;
    return FALSE;
}

void output_pbws(struct perBaseWig *pbw_list, struct slDouble *c_list, int decimals, enum wigOutType wot, boolean skip_NA, boolean skip_min, double min, FILE *out)
/* outputs one set of perBaseWigs all at the same section */
{
    struct perBaseWig *pbw;
    if (pbw_list)
    {
	boolean last_skipped = TRUE;
	int i = 0;
	int last_printed = -2;
	for (i = 0; i < pbw_list->len; i++)
	{
	    if ((!skip_NA || !has_na(pbw_list, i)) && (!skip_min || !has_under(pbw_list, i, min)))
	    {
		if (i - last_printed > 1)
		{
		    if (wot == varStepOut)
			fprintf(out, "variableStep chrom=%s span=1\n", pbw_list->chrom);
		    else if (wot == fixStepOut)
			fprintf(out, "fixedStep chrom=%s start=%d step=1 span=1\n", pbw_list->chrom, pbw_list->chromStart+i+1);
		}
		print_line(pbw_list, c_list, decimals, wot, i, out);
		last_printed = i;
	    }
	}
    }
}

struct slDouble *parse_constants(char *consts)
/* simply process the comma-list of constants from the command and return the list*/
{
    if (!consts)
	return NULL;
    struct slName *strings = slNameListFromComma(consts);
    struct slName *s;
    struct slDouble *c_list = NULL;
    for (s = strings; s != NULL; s = s->next)
    {
	struct slDouble *d = slDoubleNew(sqlDouble(s->name));
	slAddHead(&c_list, d);
    }
    slReverse(&c_list);
    slFreeList(&strings);
    return c_list;
}

void bwtool_paste(struct hash *options, char *favorites, char *regions, unsigned decimals, double fill,
		  enum wigOutType wot, struct slName **p_files, char *tmp_dir, char *output_file)
/* bwtool_paste - main for paste program */
{
    struct metaBig *mb;
    struct metaBig *mb_list = NULL;
    struct bed *bed;
    struct slName *file;
    int num_sections = 0;
    int i = 0;
    boolean skip_na = (hashFindVal(options, "skip-NA") != NULL) ? TRUE : FALSE;
    if (!isnan(fill) && skip_na)
	errAbort("cannot use -skip_na with -fill");
    boolean header = (hashFindVal(options, "header") != NULL) ? TRUE : FALSE;
    boolean verbose = (hashFindVal(options, "verbose") != NULL) ? TRUE : FALSE;
    boolean do_mean_consts = (hashFindVal(options, "consts-means") != NULL) ? TRUE : FALSE;
    boolean do_total_consts = (hashFindVal(options, "consts-totals") != NULL) ? TRUE : FALSE;
    boolean do_covs_consts = (hashFindVal(options, "consts-covs") != NULL) ? TRUE : FALSE;
    boolean skip_min = FALSE;
    double min = 0;
    if (hashFindVal(options, "skip-min"))
    {
	skip_min = TRUE;
	char *min_s = (char *)hashFindVal(options, "skip-min");
	if (min_s)
	    min = sqlDouble(min_s);
	else
	    min = 0;
    }
    struct slDouble *c_list = parse_constants((char *)hashOptionalVal(options, "consts", NULL));
    struct slDouble *fix_consts = NULL;
    struct slName *labels = NULL;
    struct slName *files = *p_files;
    FILE *out = (output_file) ? mustOpen(output_file, "w") : stdout;
    /* open the files one by one */
    if (slCount(files) == 1)
	check_for_list_files(&files, &labels, 0);
    for (file = files; file != NULL; file = file->next)
    {
	mb = metaBigOpenWithTmpDir(file->name, tmp_dir, regions);
	if (do_mean_consts || do_total_consts || do_covs_consts)
	{
	    struct bbiSummaryElement sum = bbiTotalSummary(mb->big.bbi);
	    if (do_mean_consts)
	    {
		struct slDouble *d = slDoubleNew((double)sum.sumData/sum.validCount);
		slAddHead(&fix_consts, d);
	    }
	    if (do_total_consts)
	    {
		struct slDouble *d = slDoubleNew((double)sum.sumData);
		slAddHead(&fix_consts, d);
	    }
	    if (do_covs_consts)
	    {
		struct slDouble *d = slDoubleNew((double)sum.validCount);
		slAddHead(&fix_consts, d);
	    }
	}
	slAddHead(&mb_list, mb);
    }
    slReverse(&mb_list);
    if (fix_consts)
    {
	slReverse(&fix_consts);
	c_list = slCat(c_list, fix_consts);
    }
    num_sections = slCount(mb_list->sections);
    if (header)
    {
	printf("#chrom\tchromStart\tchromEnd");
	if (labels)
	{
	    struct slName *label;
	    for (label = labels; label != NULL; label = label->next)
		printf("\t%s", label->name);
	}
	else
	    for (mb = mb_list; mb != NULL; mb = mb->next)
		printf("\t%s", mb->fileName);
	if (c_list)
	{
	    int i;
	    int size = slCount(c_list);
	    for (i = 0; i < size; i++)
		printf("\tConstant_%d", i+1);
	}
	printf("\n");
    }
    for (bed = mb_list->sections; bed != NULL; bed = bed->next)
    {
	struct perBaseWig *pbw_list = NULL;
	/* load each region */
	if (verbose)
	    fprintf(stderr, "section %d / %d: %s:%d-%d\n", i++, num_sections, bed->chrom, bed->chromStart, bed->chromEnd);
	for (mb = mb_list; mb != NULL; mb = mb->next)
	{
	    struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, bed->chrom, bed->chromStart, bed->chromEnd, FALSE, fill);
	    /* if the load returns null then NA the whole thing. */
	    /* this isn't very efficient but it's the easy way out. */
	    if (!pbw)
		pbw = alloc_perBaseWig(bed->chrom, bed->chromStart, bed->chromEnd);
	    slAddHead(&pbw_list, pbw);
	}
	slReverse(&pbw_list);
	output_pbws(pbw_list, c_list, decimals, wot, skip_na, skip_min, min, out);
	perBaseWigFreeList(&pbw_list);
    }
    /* close the files */
    carefulClose(&out);
    while ((mb = slPopHead(&mb_list)) != NULL)
	metaBigClose(&mb);
    if (labels)
	slNameFreeList(&labels);
    slNameFreeList(p_files);
    if (c_list)
	slFreeList(&c_list);
}
