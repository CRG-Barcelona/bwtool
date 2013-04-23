/* bwtool_matrix - with a bed file of regions, select  */

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

void usage_matrix()
/* Explain usage of matrix program and exit. */
{
errAbort(
  "bwtool matrix - extract data and output in a tab-delimited way that is\n"
  "   easily used as a matrix by other programs.\n"
  "usage:\n"
  "   bwtool matrix regions.bed input.bw output.txt\n"
  "options:\n"
  "   -keep-bed       in this case output the original bed loci in the first\n"
  "                   columns of the output and output data as comma-separated\n"
  "   -cluster=k      cluster regions with k-means where k is the number of\n"
  "                   clusters\n"
  "   -cluster-centroids=file\n"
  "                   store the calculated cluster centroids in a file additional\n"
  "                   to output.txt\n"
  );
}

void bwtool_matrix(struct hash *options, char *favorites, char *regions, unsigned decimals, 
		   char *bigfile, char *outputfile)
/* bwtool_matrix - main for matrix-creation program */
{
    struct metaBig *mb = metaBigOpen_favs(bigfile, regions, favorites);
    FILE *out = mustOpen(outputfile, "w");
    boolean do_k = (hashFindVal(options, "cluster") != NULL) ? TRUE : FALSE;
    boolean keep_bed = (hashFindVal(options, "keep-bed") != NULL) ? TRUE : FALSE;
    int k = (int)sqlUnsigned((char *)hashOptionalVal(options, "cluster", "0"));
    if ((do_k) && ((k < 2) || (k > 10)))
	errAbort("k should be between 2 and 10\n");
    struct bed6 *regs = readBed6(regions);
    struct perBaseMatrix *pbm = load_perBaseMatrix(mb, regs);
    struct cluster_bed_matrix *cbm = NULL;
    if (do_k)
    {
	/* ordered by cluster with label in first column */
	int i, j;
	cbm = init_cbm_from_pbm(pbm, k);
	do_kmeans(cbm, 0.001);
	for (i = 0; i < cbm->pbm->nrow; i++)
	{
	    struct perBaseWig *pbw = cbm->pbm->array[i];
	    if (keep_bed)
		fprintf(out, "%s\t%d\t%d\t%s\t%d\t%c\t", pbw->chrom, pbw->chromStart, pbw->chromEnd, pbw->name, pbw->score, pbw->strand[0]);
	    fprintf(out, "%d\t", pbw->label);
	    for (j = 0; j < pbw->len; j++)
		fprintf(out, "%0.*f%c", decimals, pbw->data[j], (j == pbw->len-1) ? '\n' : '\t');
	}
	free_cbm(&cbm);
    }
    else 
    {
	/* unordered, no label  */
	int i,j;
	for (i = 0; i < pbm->nrow; i++)
	{
	    struct perBaseWig *pbw = pbm->array[i];
	    if (keep_bed)
		fprintf(out, "%s\t%d\t%d\t%s\t%d\t%c\t", pbw->chrom, pbw->chromStart, pbw->chromEnd, pbw->name, pbw->score, pbw->strand[0]);
	    for (j = 0; j < pbw->len; j++)
		fprintf(out, "%0.*f%c", decimals, pbw->data[j], (j == pbw->len-1) ? '\n' : '\t');
	}	
	free_perBaseMatrix(&pbm);
    }
    bed6FreeList(&regs);
    carefulClose(&out);
    metaBigClose(&mb);
}
