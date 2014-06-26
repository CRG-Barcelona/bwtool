/* Some common routines for a few of the bwtool programs. */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <jkweb/common.h>
#include <jkweb/sqlNum.h>
#include <beato/metaBig.h>
#include <beato/bigs.h>
#include <jkweb/bigWig.h>
#include <jkweb/bwgInternal.h>
#include "bwtool_shared.h"

#include <math.h>

struct bed6 *load_and_recalculate_coords(char *list_file, int left, int right, boolean firstbase, boolean starts, boolean ends)
/* do the coordinate recalculation */
{
    struct bed6 *bed;
    struct bed6 *list = readBed6Soft(list_file);
    for (bed = list; bed != NULL; bed = bed->next)
    {
	boolean rev = (bed->strand[0] == '-');
	int reg_size = bed->chromEnd - bed->chromStart;
	boolean odd_sized = (reg_size % 2 == 1);
	int center = 0;
	if (!rev)
	{
	    if (starts)
		center = bed->chromStart;
	    else if (ends)
		center = bed->chromEnd;
	    else
		center = bed->chromStart + (reg_size/2);
	    bed->chromStart = center - left;
	    bed->chromEnd = center + right;
	    if (firstbase)
	    {
		if (ends)
		    bed->chromStart--;
		else
		    bed->chromEnd++;
	    }
	}
	else
	{
	    if (starts)
		center = bed->chromEnd;
	    else if (ends)
		center = bed->chromStart;
	    else if (!odd_sized)
		center = bed->chromStart + (reg_size/2);
	    else
		center = bed->chromStart + (reg_size/2) + 1;
	    bed->chromStart = center - right;
	    bed->chromEnd = center + left;
	    if (firstbase)
	    {
		if (ends)
		    bed->chromEnd++;
		else
		    bed->chromStart--;
	    }
	}
    }
    return list;
}

static boolean sniff_bed6(char *filename)
/* read up to 10 lines of the file and if nothing's wrong, return true */
{
    char *words[6];
    struct lineFile *lf = lineFileOpen(filename, TRUE);
    int i = 0;
    int num_words = 0;
    boolean good_bed = TRUE;
    while (((num_words = lineFileChop(lf, words)) == 6) && (i < 10))
    {
	if (((words[5][0] != '+') && (words[5][0] != '-')) || (countLeadingDigits(words[1]) != strlen(words[1]))
	    || (countLeadingDigits(words[2]) != strlen(words[2])))
	{
	    good_bed = FALSE;
	    break;
	}
	i++;
    }
    lineFileClose(&lf);
    if (num_words != 0)
	return FALSE;
    return good_bed;
}

static struct slName *possibly_read_list(char *filename, struct slName **pLabel_list, int ix)
/* Hopefully the file is a list of files.  return this, with labels set.  */
{
    if (isBigWigOrBed(filename) != isNotBig)
	return NULL;
    if (sniff_bed6(filename) == TRUE)
	return NULL;
    struct lineFile *lf = lineFileOpen(filename, TRUE);
    char *words[6];
    int num_words;
    struct slName *file_list = NULL;
    struct slName *label_list = NULL;
    int line_num = 1;
    while ((num_words = lineFileChop(lf, words)) > 0)
    {
	if (fileExists(words[0]))
	{
	    if (!ix || (ix == line_num))
	    {
		struct slName *newfile = slNameNew(words[0]);
		slAddHead(&file_list, newfile);
		if (num_words > 1)
		{
		    struct slName *label = slNameNew(words[1]);
		    slAddHead(&label_list, label);
		}
	    }
	    line_num++;
	}
    }
    if (line_num < ix)
	errAbort("file number %d is too high.  only %d files found.", ix, line_num);
    slReverse(&file_list);
    slReverse(&label_list);
    if (pLabel_list != NULL)
	*pLabel_list = label_list;
    return file_list;
}

int check_for_list_files(struct slName **pList, struct slName **lf_list_labels, int ix)
/* check if the comma-lists contain actual files or if they contain files listing files */
/* either way, it should expand the file into a list and return the size of the list */
/* if it's the case that the comma-list are files, then nothing should change. */
/* ix only works if it's > 0 and there isn't multiple comma-separated files */
{
    struct slName *list = *pList;
    struct slName *tmp_list = NULL;
    struct slName *label_list = NULL;
    struct slName *cur;
    int ixx = 0;
    if ((slCount(list) == 1) && (ix > 0))
	ixx = ix;
    if (!list)
	errAbort("Something went terribly wrong trying to load a file");
    while ((cur = slPopHead(&list)) != NULL)
    {
	struct slName *potential_lf_labels = NULL;
	struct slName *potential_list = possibly_read_list(cur->name, &potential_lf_labels, ixx);
	if (potential_list)
	{
	    struct slName *cat = (tmp_list) ? slCat(tmp_list, potential_list) : potential_list;
	    tmp_list = cat;
	    cat = (label_list) ? slCat(label_list, potential_lf_labels) : potential_lf_labels;
	    label_list = cat;
	}
	else
	    slAddTail(&tmp_list, cur);
    }
    if (lf_list_labels != NULL)
	*lf_list_labels = label_list;
    *pList = tmp_list;
    return slCount(tmp_list);
}

int parse_left_right(char *size_s, unsigned *pleft, unsigned *pright, int *pmeta)
/* parse the "left:right" or "left:meta:right" from the command */
/* return the number of args: 2 or 3 */
{
    unsigned left = 0, right = 0;
    int meta = 0;
    char *range[3];
    int range_num = chopString(size_s, ":", range, 3);
    if ((range_num != 2) && (range_num != 3))
	errAbort("wrongly formatted range up:down or up:meta:down");
    left = sqlUnsigned(range[0]);
    *pleft = left;
    if (range_num == 2)
    {
	right = sqlUnsigned(range[1]);
	*pright = right;
    }
    else
    {
	meta = sqlSigned(range[1]);
	*pmeta = meta;
	right = sqlUnsigned(range[2]);
	*pright = right;
    }
    return range_num;
}

void writeBw(char *inName, char *outName, struct hash *chromSizeHash)
/* shared func */
{
    struct lm *lm = lmInit(0);
    struct bwgSection *sectionList = bwgParseWig(inName, TRUE, chromSizeHash, 1024, lm);
    if (sectionList == NULL)
	errAbort("%s is empty of data", inName);
    bwgCreate(sectionList, chromSizeHash, 256, 1024, TRUE, outName);
    lmCleanup(&lm);
}

static boolean local_file(char *filename)
/* return TRUE if the file is avialable locally */
{
    char *s = cloneString(filename);
    char *colon = strchr(s, ':');
    boolean ret = FALSE;
    if (colon)
	*colon = '\0';
    if (fileExists(s))
	ret = TRUE;
    freeMem(s);
    return ret;
}

struct metaBig *metaBigOpen_check(char *bigfile, char *tmp_dir, char *regions)
/* A wrapper for metaBigOpen that does some checking and erroring */
{
    struct metaBig *mb = metaBigOpenWithTmpDir(bigfile, tmp_dir, regions);
    if (!mb)
    {
	boolean internet = (strstr(bigfile, "tp://") != 0);
	if (!internet && !local_file(bigfile))
	    errAbort("%s wasn't found", bigfile);
	else if (!internet)
	    errAbort("%s could not be opened. Perhaps it is not a valid bigWig.", bigfile);
	else
	    errAbort("There was a problem opening %s. Perhaps there is a problem with the internet connection.", bigfile);
    }
    return mb;
}

void fuse_pbm(struct perBaseMatrix **pBig, struct perBaseMatrix **pTo_add, boolean add_coords)
/* not this makes perhaps-illegal perBaseWigs where the chromEnd-chromStart are not the */
/* same as the len... which may break things somewhere if this were ever library-ized */
{
    if (pBig && pTo_add && *pTo_add)
    {
	if (*pBig == NULL)
	{
	    *pBig = *pTo_add;
	}
	else
	{
	    struct perBaseMatrix *big = *pBig;
	    struct perBaseMatrix *to_add = *pTo_add;
	    if (to_add->nrow == big->nrow)
	    {
		int i, j;
		int expanding = to_add->ncol;
		for (i = 0; i < big->nrow; i++)
		{
		    struct perBaseWig *big_pbw = big->array[i];
		    struct perBaseWig *add_pbw = to_add->array[i];
		    struct perBaseWig *new_pbw = alloc_perBaseWig(big_pbw->chrom, big_pbw->chromStart, big_pbw->chromStart + big->ncol + expanding);
		    new_pbw->name = cloneString(big_pbw->name);
		    new_pbw->score = 0;
		    new_pbw->strand[0] = big_pbw->strand[0];
		    new_pbw->chromEnd = (add_coords) ? add_pbw->chromEnd : big_pbw->chromEnd;
		    for (j = 0; j < big_pbw->len; j++)
			new_pbw->data[j] = big_pbw->data[j];
		    for (j = 0; j < add_pbw->len; j++)
			new_pbw->data[j+big_pbw->len] = add_pbw->data[j];
		    big->array[i] = new_pbw;
		    big->matrix[i] = new_pbw->data;
		    perBaseWigFree(&big_pbw);
		}
		big->ncol += expanding;
		free_perBaseMatrix(pTo_add);
	    }
	}
    }
}

int calculate_meta_file(char *file_name)
/* from all the beds in all the region files, get a single average */
{
    int count = 0;
    int sum = 0;
    if (!file_name)
	return 0;
    struct bed6 *beds = readBed6Soft(file_name);
    struct bed6 *bed;
    for (bed = beds; bed != NULL; bed = bed->next)
    {
	count++;
	sum += bed->chromEnd - bed->chromStart;
    }
    bed6FreeList(&beds);
    return sum / count;
}

int calculate_meta_file_list(struct slName *region_list)
/* from all the beds in all the region files, get a single average */
{
    int count = 0;
    int sum = 0;
    struct slName *reg;
    if (!region_list)
	return 0;
    for (reg = region_list; reg != NULL; reg = reg->next)
    {
	struct bed6 *beds = readBed6Soft(reg->name);
	struct bed6 *bed;
	for (bed = beds; bed != NULL; bed = bed->next)
	{
	    count++;
	    sum += bed->chromEnd - bed->chromStart;
	}
	bed6FreeList(&beds);
    }
    return sum / count;
}

