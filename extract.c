/* bwtool_extract - extrancting data in other ways than matrix, paste, or window */

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

#include <math.h>

#define NANUM sqrt(-1)

void usage_extract()
/* Explain usage and exit. */
{
errAbort(
  "bwtool extract - extract data from the bigWig in other ways than matrix, paste, or window\n"
  "usage:\n"
  "   bwtool extract <style> regions.bed in.bw out.txt\n"
  "where \"style\" is one of:\n"
  "   bed  - this will do something similar to bwtool matrix without the left:right specif-\n"
  "          ication and data only coming from the defined bed region, meaning region sizes\n"
  "          are also allowed to be variably-sized.  If a six-field bed is given and the\n"
  "          region is on the minus strand, then the extracted data is reversed prior to\n"
  "          outputting. The output format is the original bed up to the first six fields,\n"
  "          tab-delimited, followed by a field indicating the length of the data to follow,\n"
  "          followed by the data, separated by commas (or tabs if option -tabs is used).\n"
  "   jsp  - with similar effect as \"bed\" in terms of stranded bed input and reversing data\n"
  "          or not, the output is a bit more minimal and has a vertical structure simlar to\n"
  "          bwtool paste.  Values are separated line-by-line and regions are preceded by a\n"
  "          line starting with # and stating the name of the region from the bed. If only\n"
  "          three fields are used in the bed or bed name is \".\", then the region is numbered.\n"
  "options:\n"
  "   -tabs         output tabs instead of commas in output.\n"
  "   -locus-name   in jsp output, output the region locus in genome browser coordinate\n"
  "                 form (i.e. chrom:(chromStart+1)-chromEnd instead of the bed name\n"
  );
}

enum style_type
{
    nothing = 0,
    bed = 1,
    jsp = 2,
};

void extractOutBed(FILE *out, struct bed6 *section, int orig_size, unsigned decimals, struct perBaseWig *pbw, boolean tabs)
/* Do the output like:
   chr1   2   4   2   3.00,4.00
   chr1   8   12  4   9.00,10.00,11.00,12.00    */
{
    int i;
    fprintf(out, "%s\t%d\t%d\t", section->chrom, section->chromStart, section->chromEnd);
    if (orig_size > 3)
	fprintf(out, "%s\t", section->name);
    if (orig_size > 4)
	fprintf(out, "%d\t", section->score);
    if (orig_size > 5)
	fprintf(out, "%c\t", section->strand[0]);
    fprintf(out, "%d\t", pbw->len);
    for (i = 0; i < pbw->len-1; i++)
	if (isnan(pbw->data[i]))
	    fprintf(out, "NA%c", (tabs) ? '\t' : ',');
	else
	    fprintf(out, "%0.*f%c", decimals, pbw->data[i], (tabs) ? '\t' : ',');
    if (isnan(pbw->data[pbw->len-1]))
	fprintf(out, "NA\n");
    else
	fprintf(out, "%0.*f\n", decimals, pbw->data[pbw->len-1]);
}

void extractOutJsp(FILE *out, struct bed6 *section, unsigned decimals, struct perBaseWig *pbw)
/* Do the output like:
   # region_1
   3.00
   4.00
   # region_2
   12.00
   11.00
   10.00
   9.00              */
{
    int i;
    fprintf(out, "# %s\n", section->name);
    for (i = 0; i < pbw->len; i++)
	if (isnan(pbw->data[i]))
	    fprintf(out, "NA\n");
	else
	    fprintf(out, "%0.*f\n", decimals, pbw->data[i]);
}

void bwtool_extract(struct hash *options, char *regions, unsigned decimals, double fill,
		  char *style_s, char *bigfile, char *tmp_dir, char *outputfile)
/* bwtool_extract - main for the extract program */
{
    boolean tabs = (hashFindVal(options, "tabs") != NULL) ? TRUE : FALSE;
    boolean locus_name = (hashFindVal(options, "locus-name") != NULL) ? TRUE : FALSE;
    int orig_size = 0;
    struct bed6 *region_list = readBed6SoftAndSize(regions, &orig_size);
    struct metaBig *mb = metaBigOpenWithTmpDir(bigfile, tmp_dir, NULL);
    if (!mb)
	errAbort("problem opening %s", bigfile);
    FILE *out = mustOpen(outputfile, "w");
    struct bed6 *section;
    enum style_type style = nothing;
    if (sameWord(style_s, "bed"))
	style = bed;
    else if (sameWord(style_s, "jsp"))
	style = jsp;
    else
	errAbort("please specify a valid style");
    int section_num = 1;
    /* loop through each region */
    for (section = region_list; section != NULL; section = section->next)
    {
	struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, section->chrom, section->chromStart,
							      section->chromEnd, (section->strand[0] == '-') ? TRUE : FALSE, fill);
	if (style == bed)
	    /* for bed there is no name manipulation */
	    extractOutBed(out, section, orig_size, decimals, pbw, tabs);
	else
	{
	    /* for jsp output there is some name manipulation that could be done prior to outputting */
	    char buf[128];
	    if (locus_name)
	    {
		safef(buf, sizeof(buf), "%s:%d-%d", section->chrom, section->chromStart+1, section->chromEnd);
		if (section->name)
		    freeMem(section->name);
		section->name = cloneString(buf);
	    }
	    else if ((orig_size < 4) || sameString(section->name, "."))
	    {
		safef(buf, sizeof(buf), "region_%d", section_num);
		if (section->name)
		    freeMem(section->name);
		section->name = cloneString(buf);
	    }
	    extractOutJsp(out, section, decimals, pbw);
	}
	perBaseWigFree(&pbw);
	section_num++;
    }
    metaBigClose(&mb);
    bed6FreeList(&region_list);
    carefulClose(&out);
}
