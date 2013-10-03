#include "common.h"
#include "metaBig.h"
#include "bigs.h"
#include "bwtool_shared.h"

struct bed6 *load_and_recalculate_coords(char *list_file, int left, int right, boolean firstbase, boolean starts, boolean ends)
/* do the coordinate recalculation */
{
    struct bed6 *bed;
    struct bed6 *list = readBed6Soft(list_file);
    for (bed = list; bed != NULL; bed = bed->next)
    {
	boolean rev = (bed->strand[0] == '-');
	int l = (rev) ? right : left;
	int size = left + right + (firstbase ? 1 : 0);
	int center = bed->chromStart + ((bed->chromEnd - bed->chromStart) / 2);
	if (starts || ends)
	{
	if ((starts && !rev) || (ends && rev))
	    center = bed->chromStart;
	else if ((ends && !rev) || (starts && rev))
	    center = bed->chromEnd;
	}
	bed->chromStart = center - l;
	bed->chromEnd = bed->chromStart + size;
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
