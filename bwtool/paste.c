/* bwtool_paste - simultaneously output same regions of multiple bigWigs */

#include "common.h"
#include "obscure.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "basicBed.h"
#include "bigWig.h"
#include "bigs.h"
#include "bwtool.h"
#include "cluster.h"

void usage_paste()
/* Explain usage of paste program and exit. */
{
errAbort(
  "bwtool paste - simultaneously output same regions of multiple bigWigs\n"
  "usage:\n"
  "   bwtool paste input1.bw input2.bw input3.bw ...\n"
  "options:\n"
  "   -skip-NA        don't output lines (bases) where one of the inputs is NA\n"
  );
}

/* void check_chrom_sizes(struct metaBig *mbList) */
/* /\* check to make sure the metaBigs all have the same size chroms. *\/ */
/* /\* After all, it would be bad to mix genomes/assemblies by accident. *\/ */
/* { */
/*     /\* First trim down all chroms to one list *\/ */
/*     if (!mbList) */
/* 	errAbort("No metaBigs in input"); */
/*     struct hash *one_hash = mbList->chromSizeHash; */
/*     struct metaBig *mb; */
/*     struct hashEl *list = hashElListHash(one_hash); */
/*     for (mb = mbList; mb != NULL; mb = mb->next) */
/*     { */
/* 	struct hashEl *el; */
/* 	for (el =  */
/*     } */
/* } */

void output_pbws(struct perBaseWig *pbw_list, int decimals, boolean skip_NA)
/* outputs one set of perBaseWigs all at the same section */
{
    struct perBaseWig *pbw;
    if (pbw_list)
    {
	int i = 0;
	while (i < pbw_list->len)
	{
	    if (skip_NA)
	    {
		for (pbw = pbw_list; pbw != NULL; pbw = pbw->next)
		    if (isnan(pbw->data[i]))
		    {
			i++;
			continue;
		    }
	    }
	    printf("%s\t%d\t%d\t", pbw_list->chrom, pbw_list->chromStart+i, pbw_list->chromStart+i+1);
	    for (pbw = pbw_list; pbw != NULL; pbw = pbw->next)
	    {
		if (isnan(pbw->data[i]))
		    printf("NA");
		else
		    printf("%0.*f", decimals, pbw->data[i]);
		printf("%c", (pbw->next == NULL) ? '\n' : '\t');
	    }
	    i++;
	}
    }
}

void bwtool_paste(struct hash *options, char *favorites, char *regions, unsigned decimals, 
		   struct slName *files)
/* bwtool_paste - main for paste program */
{
    struct metaBig *mb;
    struct metaBig *mb_list = NULL;
    struct bed6 *reg_list = readBed6(regions);
    struct bed6 *bed; 
    struct slName *file;
    boolean skip_na = (hashFindVal(options, "skip-NA") != NULL) ? TRUE : FALSE;
    /* open the files one by one */
    for (file = files; file != NULL; file = file->next)
    {
	mb = metaBigOpen(file->name, regions);
	slAddHead(&mb_list, mb);
    }
    /* list is reversed but then so will be making the list of pbws, */
    /* so this avoids double-reversing */
    for (bed = mb->sections; bed != NULL; bed = bed->next)
    {
	struct perBaseWig *pbw_list = NULL;
	/* load each region */
	for (mb = mb_list; mb != NULL; mb = mb->next)
	{
	    struct perBaseWig *pbw = perBaseWigLoadSingleContinue(mb, bed->chrom, bed->chromStart, bed->chromEnd, FALSE);
	    slAddHead(&pbw_list, pbw);
	}
	output_pbws(pbw_list, decimals, skip_na);
	perBaseWigFreeList(&pbw_list);
    }
    /* close the files */
    while ((mb = slPopHead(&mb_list)) != NULL)
	metaBigClose(&mb);
}
