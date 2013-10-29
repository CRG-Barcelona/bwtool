#include "common.h"
#include "sqlNum.h"
#include "metaBig.h"
#include "bigs.h"
#include "bigWig.h"
#include "bwgInternal.h"
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

static struct slName *possibly_read_list(char *filename, struct slName **pLabel_list)
/* */
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
    while ((num_words = lineFileChop(lf, words)) > 0)
    {
	if (fileExists(words[0]))
	{
	    struct slName *newfile = slNameNew(words[0]);
	    slAddHead(&file_list, newfile);
	    if (num_words > 1)
	    {
		struct slName *label = slNameNew(words[1]);
		slAddHead(&label_list, label);
	    }
	}
    }
    slReverse(&file_list);
    slReverse(&label_list);
    if (pLabel_list != NULL)
	*pLabel_list = label_list;
    return file_list;
}

int check_for_list_files(struct slName **pList, struct slName **lf_list_labels)
/* check if the comma-lists contain actual files or if they contain files listing files */
/* either way, it should expand the file into a list and return the size of the list */
/* if it's the case that the comma-list are files, then nothing should change. */
{
    struct slName *list = *pList;
    struct slName *tmp_list = NULL;
    struct slName *label_list = NULL;
    struct slName *cur;
    if (!list)
	errAbort("Something went terribly wrong trying to load a file");
    while ((cur = slPopHead(&list)) != NULL)
    {
	struct slName *potential_lf_labels = NULL;
	struct slName *potential_list = possibly_read_list(cur->name, &potential_lf_labels);
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

void parse_left_right(char *size_s, unsigned *pleft, unsigned *pright)
/* parse the "left:right" from the command */
{
    unsigned left = 0, right = 0;
    char *tmp_s = cloneString(size_s);
    char *range[2];
    int range_num = chopString(tmp_s, ":", range, sizeof(range));
    if (range_num != 2)
	errAbort("wrongly formatted range left:right");
    left = sqlUnsigned(range[0]);
    right = sqlUnsigned(range[1]);
    *pleft = left;
    *pright = right;
    freeMem(tmp_s);
}

void writeBw(char *inName, char *outName, struct hash *chromSizeHash)
/* shared func */
{
    struct lm *lm = lmInit(0);
    struct bwgSection *sectionList = bwgParseWig(inName, TRUE, chromSizeHash, 1024, lm);
    if (sectionList == NULL)
	errAbort("%s is empty of data", inName);
    bwgCreate(sectionList, chromSizeHash, 256, 1024, FALSE, outName);
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

struct metaBig *metaBigOpen_check(char *bigfile, char *regions)
/* A wrapper for metaBigOpen that does some checking and erroring */
{
    struct metaBig *mb = metaBigOpen(bigfile, regions);
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
