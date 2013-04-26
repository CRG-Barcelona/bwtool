/* bwtool_centrplot - make the typical centerplots */

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "basicBed.h"
#include "bigWig.h"
#include "bigs.h"
#include "random_coord.h"
#include "bwtool.h"
#include "cluster.h"
#include "stuff.h"

#define NANUM sqrt(-1)

void usage_aggregate()
/* Explain usage of distribution program and exit. */
{
errAbort(
  "bwtool aggregate - produce plot data as averages surrounding given regions\n"
  "   seen in the bigWig.\n"
  "usage:\n"
  "   bwtool aggregate left:right a.bed,b.bed,... x.bw,y.bw,... output.txt\n" 
  "where:\n"
  "  left:right is a description of the range to use surrounding the\n"
  "         centers/starts. A value of 200:300 will create plot points\n"
  "         for -200 bp to +300 bp from the center of all the regions.\n"
  "  a comma-separated list of bed files (or just one) OR a comma-\n"
  "         separated list of data bigWigs (or just one) is given\n"
  "options:\n"
  "   -starts          use starts of bed regions as opposed to the middles\n"
  "   -ends            use ends of bed regions as opposed to the middles\n"
  "   -nozero          skip zero in the output i.e. -3,-2,-1,1,2,3,\n"
  "                    this has the effect of the limits being\n"
  "                    -size to size instead of -size to (size-1)\n"
  "   -firstbase       in this case the zero base is used so output\n"
  "                    has 2*size+1 lines (incompatible with -nozero)\n"
  "   -expanded        output medians and standard deviations instead of just\n"
  "                    averages\n" 
  "   -cluster=k       cluster with k-means with given parameter (2-10 are best)\n"
  "   -cluster-sets=file.bed\n"
  "                    write out the original bed with the cluster label along\n"
  "                    with the accompanying region used for the clustering\n"
  "   -verbose=n       output some progress info every n regions loaded\n"
  );
}

struct agg_data
{
    int nrow;
    int ncol;
    int *indexes;
    double **data;
};

struct agg_data *init_agg_data(int left, int right, boolean firstbase, boolean nozero, int num_cols)
/* init the output struct */
{
    struct agg_data *agg;
    int i;
    int addone = 0;
    AllocVar(agg);
    agg->nrow = left + right;
    agg->ncol = num_cols;
    if (firstbase)
	agg->nrow++;
    AllocArray(agg->indexes, agg->nrow);
    for (i = 0; i < left; i++)
	agg->indexes[i] = i - left;
    if (nozero)
	addone = 1;
    for (; i < agg->nrow; i++)
	agg->indexes[i] = i - left + addone;
    AllocArray(agg->data, agg->nrow);
    for (i = 0; i < agg->nrow; i++)
	AllocArray(agg->data[i], agg->ncol);
    return agg;
}

void free_agg_data(struct agg_data **pAgg)
/* free the output struct */
{
    struct agg_data *agg = *pAgg;
    int i;
    for (i = 0; i < agg->nrow; i++)
	freeMem(agg->data[i]);
    freeMem(agg->data);
    freeMem(agg->indexes);
    freez(&agg);
}

void do_summary(struct perBaseMatrix *pbm, struct agg_data *agg, boolean expanded, int offset)
/* calculate mean, median, sd */
{
    const double na = NANUM;
    int i, j;
    double *one_pbm_col; 
    AllocArray(one_pbm_col, pbm->nrow);
    for (i = 0; i < agg->nrow; i++)
    {
	int size = 0;
	double sum = 0;
	double mean = 0;
	double sd = 0;
	for (j = 0; j < pbm->nrow; j++)
	    if (!isnan(pbm->matrix[j][i]))
		one_pbm_col[size++] = pbm->matrix[j][i];
	if (size > 0)
	{
	    for (j = 0; j < size; j++)
		sum += one_pbm_col[j];
	    mean = sum/size;
	    sum = 0;
	    for (j = 0; j < size; j++)
		sum += pow(one_pbm_col[j] - size,2);
	    sd = sqrt(sum/size);
	    agg->data[i][offset] = mean;
	    if (expanded)
	    {
		agg->data[i][offset+1] = doubleMedian(size, one_pbm_col);
		agg->data[i][offset+2] = sd;
		agg->data[i][offset+3] = (double)size;
	    }
	}
	else
	{
	    agg->data[i][offset] = na;
	    if (expanded)
	    {
		agg->data[i][offset+1] = na;
		agg->data[i][offset+2] = na;
		agg->data[i][offset+3] = na;
	    }
	}
    }
    freeMem(one_pbm_col);
}

void copy_centroids(struct cluster_bed_matrix *cbm, struct agg_data *agg)
/* copies the centroid values into the output struct */
{
    int i, j;
    for (i = 0; i < cbm->k; i++)
    {
	for (j = 0; j < cbm->m; j++)
	    agg->data[j][i] = cbm->centroids[i][j];
    }
}

void output_agg_data(FILE *out, int decimals, boolean expanded, struct agg_data *agg)
/* simply output the stuff */
{
    int i, j;
    for (i = 0; i < agg->nrow; i++)
    {
	fprintf(out, "%d\t", agg->indexes[i]);
	j = 0;
	while (j < agg->ncol)
	{
	    fprintf(out, "%0.*f", decimals, agg->data[i][j]);
	    if (expanded)
	    {
		fprintf(out, "\t%0.*f\t%0.*f\t%d", decimals, agg->data[i][j+1], decimals, agg->data[i][j+2], (int)agg->data[i][j+3]);
		j += 4;
	    }
	    else
		j++;
	    fprintf(out, "%c", (j < agg->ncol) ? '\t' : '\n');
	}
    }
}

void output_cluster_sets(struct cluster_bed_matrix *cbm, char *cluster_sets)
/* Output original bed6, followed by the cluster label, followed by the modified bed6 */
/* indicating the region used in clustering */
{
    FILE *out = mustOpen(cluster_sets, "w");
    struct perBaseMatrix *pbm = cbm->pbm;
    int i;
    for (i = 0; i < pbm->nrow; i++)
    {
	struct perBaseWig *pbw = pbm->array[i];
	struct bed6 *ob = pbw->orig_bed;
	fprintf(out, "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%s\t%d\t%d\t%s\t%d\t%c\n", ob->chrom, 
		ob->chromStart, ob->chromEnd, ob->name, ob->score, ob->strand[0], 
		pbw->label, pbw->chrom, pbw->chromStart, pbw->chromEnd, pbw->name, pbw->score,
		pbw->strand[0]);
    }
    carefulClose(&out);
}

struct bed6 *load_and_recalculate(char *list_file, int left, int right, boolean firstbase, boolean starts, boolean ends)
/* do the coordinate recalculation */
{
    struct bed6 *bed;
    struct bed6 *list = readBed6(list_file);
    for (bed = list; bed != NULL; bed = bed->next)
    {
	boolean rev = (bed->strand[0] == '-');
	int l = (rev) ? right : left;
	int size = left + right + (firstbase ? 1 : 0);
	int center = bed->chromStart + ((bed->chromEnd - bed->chromEnd) / 2);
	if ((starts && !rev) || (ends && rev))
	    center = bed->chromStart;
	else if ((ends && !rev) || (starts && rev))
	    center = bed->chromEnd;
	bed->chromStart = center - l;
	bed->chromEnd = bed->chromStart + size;
    }
    return list;
}

void bwtool_aggregate(struct hash *options, char *regions, unsigned decimals, 
		       char *size_s, char *region_list_s, char *wig, char *output_file)
/* aggregate - main */
{
    unsigned left = 0, right = 0;
    struct slName *region_list = slNameListFromComma(region_list_s);
    struct slName *wig_list = slNameListFromComma(wig); 
    boolean firstbase = (hashFindVal(options, "firstbase") != NULL) ? TRUE : FALSE;
    boolean nozero = (hashFindVal(options, "nozero") != NULL) ? TRUE : FALSE;
    boolean use_start = (hashFindVal(options, "starts") != NULL) ? TRUE : FALSE;
    boolean use_end = (hashFindVal(options, "ends") != NULL) ? TRUE : FALSE;
    boolean expanded = (hashFindVal(options, "expanded") != NULL) ? TRUE : FALSE;
    boolean clustering = (hashFindVal(options, "cluster") != NULL) ? TRUE : FALSE;
    char *cluster_sets = (char *)hashFindVal(options, "cluster-sets");
    int k = (int)sqlUnsigned((char *)hashOptionalVal(options, "cluster", "0"));
    /* int verbose = sqlSigned((char *)hashOptionalVal(options, "verbose", "0")); */
    FILE *output;
    char *tmp_s = cloneString(size_s);
    char *range[2];
    boolean mult_regions = (slCount(region_list) > 1);
    boolean mult_wigs = (slCount(wig_list) > 1);
    int range_num = chopString(tmp_s, ":", range, sizeof(range));
    if (use_start && use_end)
	errAbort("cannot specify both -starts and -ends");
    if (range_num != 2)
	errAbort("wrongly formatted range left:right");
    left = sqlUnsigned(range[0]);
    right = sqlUnsigned(range[1]);
    if ((clustering) && ((k < 2) || (k > 10)))
	errAbort("k should be between 2 and 10\n");
    if (mult_regions && mult_wigs)
	errAbort("can't specify more than one region set AND more than one dataset: must choose either multiple region sets or multiple datasets");  
    if ((mult_regions || mult_wigs) && clustering)
	errAbort("with clustering just specify one region list and one bigWig");
    if (firstbase && nozero)
	errAbort("-firsbase and -nozero cannot be used together");
    output = mustOpen(output_file, "w");
    if (mult_regions)
    {
	int num_regions = slCount(region_list);
	struct agg_data *agg = init_agg_data(left, right, firstbase, nozero, (expanded) ? (num_regions * 4) : num_regions);
	struct slName *reg;
	struct metaBig *mb = metaBigOpen(wig_list->name, NULL);
	int offset = 0;
	for (reg = region_list; reg != NULL; reg = reg->next)
	{	
	    struct bed6 *regions = load_and_recalculate(reg->name, left, right, firstbase, use_start, use_end);
	    struct perBaseMatrix *pbm = load_perBaseMatrix(mb, regions);
	    do_summary(pbm, agg, expanded, offset++); 
	    free_perBaseMatrix(&pbm);
	    bed6FreeList(&regions);
	}
	metaBigClose(&mb);
    }
    else if (mult_wigs)
    {
	int num_wigs = slCount(wig_list);
	struct agg_data *agg = init_agg_data(left, right, firstbase, nozero, (expanded) ? (num_wigs * 4) : num_wigs);
	struct slName *wig_name; 
	struct bed6 *regions = load_and_recalculate(region_list->name, left, right, firstbase, use_start, use_end);
	int offset = 0;
	for (wig_name = wig_list; wig_name != NULL; wig_name = wig_name->next)
	{
	    struct metaBig *mb = metaBigOpen(wig_name->name, NULL);
	    struct perBaseMatrix *pbm = load_perBaseMatrix(mb, regions);
	    do_summary(pbm, agg, expanded, offset++);
	    output_agg_data(output, decimals, expanded, agg);
	    free_perBaseMatrix(&pbm);
	    metaBigClose(&mb);
	}
	bed6FreeList(&regions);
    }
    else if (clustering)
    {
	struct agg_data *agg = init_agg_data(left, right, firstbase, nozero, k);
	struct metaBig *mb = metaBigOpen(wig_list->name, NULL);
	struct bed6 *regions = load_and_recalculate(region_list->name, left, right, firstbase, use_start, use_end);
	struct bed6 *orig_regions = readBed6(region_list->name);
	struct perBaseMatrix *pbm = load_perBaseMatrix(mb, regions);
	if (cluster_sets)
	    perBaseMatrixAddOrigRegions(pbm, orig_regions);
	struct cluster_bed_matrix *cbm = init_cbm_from_pbm(pbm, k);
	do_kmeans(cbm, k);
	copy_centroids(cbm, agg);
	output_agg_data(output, decimals, expanded, agg);
	if (cluster_sets)
	    output_cluster_sets(cbm, cluster_sets);
	free_cbm(&cbm);
	metaBigClose(&mb);
	bed6FreeList(&regions);
	free_agg_data(&agg);	
    }
    else
    /* easiest case*/
    {
	struct agg_data *agg = init_agg_data(left, right, firstbase, nozero, (expanded) ? 3 : 1);
	struct metaBig *mb = metaBigOpen(wig_list->name, NULL);
	struct bed6 *regions = load_and_recalculate(region_list->name, left, right, firstbase, use_start, use_end);
	struct perBaseMatrix *pbm = load_perBaseMatrix(mb, regions);
	do_summary(pbm, agg, expanded, 0);
	output_agg_data(output, decimals, expanded, agg);
	free_perBaseMatrix(&pbm);
	metaBigClose(&mb);
	bed6FreeList(&regions);
	free_agg_data(&agg);
    }
    carefulClose(&output);
    slNameFreeList(&region_list);
    slNameFreeList(&wig_list);
    freeMem(tmp_s);
}
