/* bwtool_matrix - with a bed file of regions, select  */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif 

#include "common.h"
#include "obscure.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "basicBed.h"
#include "bigWig.h"
#include "bigs.h"
#include "metaBig.h"
#include "bwtool.h"
#include "cluster.h"
#include "bwtool_shared.h"

void usage_matrix()
/* Explain usage of matrix program and exit. */
{
errAbort(
  "bwtool matrix - extract data and output in a tab-delimited way that is\n"
  "   easily used as a matrix by other programs.\n"
  "usage:\n"
  "   bwtool matrix left:right regions.bed input1.bw,input2... output.txt\n"
  "   bwtool matrix left:right regions.bed bigWigs.lst output.txt\n\n"
  "If more than one bigWig is specified, then the matrix from the second\n"
  "matrix is fused to the first, and the third to the second, etc. in the\n"
  "same left-to-right order as the bigWigs\n\n"
  "options:\n"
  "   -keep-bed       in this case output the original bed loci in the first\n"
  "                   columns of the output and output data as comma-separated\n"
  "   -starts         use starts of bed regions as opposed to the middles\n"
  "   -ends           use ends of bed regions as opposed to the middles\n"
  "   -tiled-averages=n\n"
  "                   break the left:right sized region into regions of n bases\n"
  "                   and average over those subregions as the matrix.\n"
  "   -long-form=labels\n"
  "   -long-form-header\n"
  "   -cluster=k      cluster regions with k-means where k is the number of\n"
  "                   clusters\n"
  "   -cluster-centroids=file\n"
  "                   store the calculated cluster centroids in a file additional\n"
  "                   to output.txt\n"
  );
}

void output_centroids(struct cluster_bed_matrix *cbm, char *centroid_file, int decimals)
{
    FILE *out2 = mustOpen(centroid_file, "w");
    int i, j, k = cbm->k;
    int total = cbm->n - cbm->num_na;
    fprintf(out2, "# num NA = %d, num total = %d\n", cbm->num_na, total);
    for (i = 0; i < k; i++)
    {
	fprintf(out2, "cluster %d\tsize = %d (%0.2f)\tvalues = ", i, cbm->cluster_sizes[i], (double)cbm->cluster_sizes[i]/total);
	for (j = 0; j < cbm->m; j++)
	    fprintf(out2, "%f%c", cbm->centroids[i][j], (j == cbm->m-1) ? '\n' : ',');
    }
    carefulClose(&out2);
}

void output_cluster_matrix(struct cluster_bed_matrix *cbm, int decimals, boolean keep_bed, char *outputfile)
{
    FILE *out = mustOpen(outputfile, "w");
    int i, j;
    for (i = 0; i < cbm->pbm->nrow; i++)
    {
	struct perBaseWig *pbw = cbm->pbm->array[i];
	if (keep_bed)
	    fprintf(out, "%s\t%d\t%d\t%s\t%d\t%c\t", pbw->chrom, pbw->chromStart, pbw->chromEnd, pbw->name, pbw->score, pbw->strand[0]);
	fprintf(out, "%d\t", pbw->label);
	fprintf(out, "%f\t", pbw->cent_distance);
	for (j = 0; j < pbw->len; j++)
	    fprintf(out, "%0.*f%c", decimals, pbw->data[j], (j == pbw->len-1) ? '\n' : '\t');
    }
    carefulClose(&out);
}

void output_cluster_matrix_long(struct cluster_bed_matrix *cbm, struct slName *labels, boolean keep_bed, char *outputfile, boolean header)
{
    FILE *out = mustOpen(outputfile, "w");
    int i, j, l;
    char **labels_array;
    int *subpos_array;
    int *cluster_row_array;
    int num_labels = slCount(labels);
    int num_subpos = cbm->pbm->ncol / num_labels;
    struct slName *label;
    AllocArray(labels_array, cbm->pbm->ncol);
    AllocArray(subpos_array, cbm->pbm->ncol);
    i = 0;
    for (label = labels; label != NULL; label = label->next)
	for (j = 0; j < num_subpos; j++)
	{
	    labels_array[i] = label->name;
	    subpos_array[i] = j+1;
	    i++;
	}
    AllocArray(cluster_row_array, cbm->pbm->nrow - cbm->num_na);
    i = 0;
    for (j = 0; j < cbm->k; j++)
	for (l = 0; l < cbm->cluster_sizes[j]; l++)
	    cluster_row_array[i++] = l+1;
    if (header)
    {
	if (keep_bed)
	    fprintf(out, "chrom\tchromStart\tchromEnd\tname\tscore\tstrand\t");
	fprintf(out, "Row\tCluster_Name\tCluster_Row\tCentroid_Distance\tLabel_Subpos\tLabel\tSubpos\tColumn\tData\n");
    }
    for (i = cbm->num_na; i < cbm->pbm->nrow; i++)
    {
	struct perBaseWig *pbw = cbm->pbm->array[i];
	for (j = 0; j < cbm->pbm->ncol; j++)
	{
	    if (keep_bed)
		fprintf(out, "%s\t%d\t%d\t%s\t%d\t%c\t", pbw->chrom, pbw->chromStart, pbw->chromEnd, pbw->name, pbw->score, pbw->strand[0]);
	    fprintf(out, "%d\t", i - cbm->num_na + 1);
	    fprintf(out, "Cluster_%d\t", pbw->label + 1);
	    fprintf(out, "%d\t", cluster_row_array[i-cbm->num_na]);
	    fprintf(out, "%f\t", pbw->cent_distance);
	    fprintf(out, "%s_%d\t%s\t%d\t", labels_array[j], subpos_array[j], labels_array[j], subpos_array[j]);
	    fprintf(out, "%d\t", j+1);
	    fprintf(out, "%f\n", cbm->pbm->matrix[i][j]);
	}
    }
    freeMem(labels_array);
    freeMem(subpos_array);
    carefulClose(&out);
}

void output_matrix(struct perBaseMatrix *pbm, int decimals, boolean keep_bed, char *outputfile)
{
    FILE *out = mustOpen(outputfile, "w");
    int i,j;
    for (i = 0; i < pbm->nrow; i++)
    {
	struct perBaseWig *pbw = pbm->array[i];
	if (keep_bed)
	    fprintf(out, "%s\t%d\t%d\t%s\t%d\t%c\t", pbw->chrom, pbw->chromStart, pbw->chromEnd, pbw->name, pbw->score, pbw->strand[0]);
	for (j = 0; j < pbw->len; j++)
	    fprintf(out, "%0.*f%c", decimals, pbw->data[j], (j == pbw->len-1) ? '\n' : '\t');
    }	
    carefulClose(&out);
}

static struct slName *setup_labels(char *long_form, struct slName *bw_list, struct slName **labels_from_bw_list)
{
    struct slName *lf_labels = NULL;
    if (long_form && !sameString(long_form, "on"))
    /* first check labels from CL */
    {
	lf_labels = slNameListFromComma(long_form);
	int num_labels = slCount(lf_labels);
	if (num_labels != slCount(bw_list))
	    errAbort("number of labels provided should equal the number of bigWigs given");
    }
    /* from file */
    else
    {
	if (*labels_from_bw_list)
	{
	    if (slCount(*labels_from_bw_list) == slCount(bw_list))
		lf_labels = *labels_from_bw_list;
	    else
		slNameFreeList(labels_from_bw_list);
	}
	if (!lf_labels)
	    lf_labels = slNameCloneList(bw_list);
    }
    return lf_labels;
}


void bwtool_matrix(struct hash *options, char *favorites, char *regions, unsigned decimals, 
		   double fill, char *range_s, char *bigfile, char *outputfile)
/* bwtool_matrix - main for matrix-creation program */
{
    boolean do_k = (hashFindVal(options, "cluster") != NULL) ? TRUE : FALSE;
    boolean do_tile = (hashFindVal(options, "tiled-averages") != NULL) ? TRUE : FALSE;
    boolean keep_bed = (hashFindVal(options, "keep-bed") != NULL) ? TRUE : FALSE;
    boolean starts = (hashFindVal(options, "starts") != NULL) ? TRUE : FALSE;
    boolean ends = (hashFindVal(options, "ends") != NULL) ? TRUE : FALSE;
    boolean lf_header = (hashFindVal(options, "long-form-header") != NULL) ? TRUE : FALSE;
    char *centroid_file = (char *)hashFindVal(options, "cluster-centroids");
    char *long_form = (char *)hashFindVal(options, "long-form");
    boolean do_long_form = (long_form != NULL);
    unsigned left = 0, right = 0;
    parse_left_right(range_s, &left, &right);
    int k = (int)sqlUnsigned((char *)hashOptionalVal(options, "cluster", "0"));
    int tile = (int)sqlUnsigned((char *)hashOptionalVal(options, "tiled-averages", "0"));
    if ((do_k) && ((k < 2) || (k > 10)))
	errAbort("k should be between 2 and 10\n");
    if ((do_tile) && (tile < 1))
	errAbort("tiling should be done for larger regions");
    struct slName *bw_names = slNameListFromComma(bigfile);
    struct slName *bw_name;
    struct slName *labels_from_file = NULL;
    int num_bigwigs = check_for_list_files(&bw_names, &labels_from_file);
    struct slName *labels = setup_labels(long_form, bw_names, &labels_from_file);
    struct bed6 *regs = NULL;
    struct perBaseMatrix *pbm = NULL;
    int i;
    regs = load_and_recalculate_coords(regions, left, right, FALSE, starts, ends);
    for (bw_name = bw_names; bw_name != NULL; bw_name = bw_name->next)
    {
	struct metaBig *mb = metaBigOpen(bw_name->name, NULL);
	struct perBaseMatrix *one_pbm = (do_tile) ? load_ave_perBaseMatrix(mb, regs, tile, fill) : load_perBaseMatrix(mb, regs, fill);
	fuse_pbm(&pbm, &one_pbm);
	metaBigClose(&mb);
    }
    if (do_k)
    {
	struct cluster_bed_matrix *cbm = NULL;
	/* ordered by cluster with label in first column */
	cbm = init_cbm_from_pbm(pbm, k);
	do_kmeans_sort(cbm, 0.001, TRUE);
	if (do_long_form)
	    output_cluster_matrix_long(cbm, labels, keep_bed, outputfile, lf_header);
	else
	    output_cluster_matrix(cbm, decimals, keep_bed, outputfile);
	if (centroid_file)
	    output_centroids(cbm, centroid_file, decimals);
	free_cbm(&cbm);
    }
    else 
    {
	output_matrix(pbm, decimals, keep_bed, outputfile);
	/* unordered, no label  */
	free_perBaseMatrix(&pbm);
    }
    bed6FreeList(&regs);
}
