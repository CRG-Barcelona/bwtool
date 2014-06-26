/* bwtool_centrplot - make the typical centerplots */

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
#include <beato/random_coord.h>
#include "bwtool.h"
#include "bwtool_shared.h"
#include <beato/cluster.h>
#include <beato/stuff.h>

#define NANUM sqrt(-1)
#define NUM_EXPANDED 5

void usage_aggregate()
/* Explain usage of distribution program and exit. */
{
errAbort(
  "bwtool aggregate - produce plot data as averages surrounding given regions\n"
  "   seen in the bigWig.\n"
  "usage:\n"
  "   bwtool aggregate up:down a.bed,b.bed,... x.bw,y.bw,... output.txt\n"
  "   bwtool aggregate up:down bed.lst bigWig.lst output.txt\n"
  "   bwtool aggregate up:meta:down a.bed,b.bed,... x.bw,y.bw,... output.txt\n"
  "   bwtool aggregate up:meta:down bed.lst bigWig.lst output.txt\n"
  "where:\n"
  "  up:down is a description of the range to use surrounding the\n"
  "         centers/starts. A value of 200:300 will create plot points\n"
  "         for 200 bp upstream to 300 bp downstream from the center of\n"
  "         all the regions.\n"
  "  bed.lst/bigWig.lst are files containing lists of files to be used\n"
  "         in place of a comma-separated list.\n\n"
  "options:\n"
  "   -starts          use starts of bed regions as opposed to the middles\n"
  "                    (unavailable when using meta)\n"
  "   -ends            use ends of bed regions as opposed to the middles\n"
  "                    (unavailable when using meta)\n"
  "   -firstbase       in this case the zero base is used so output\n"
  "                    has left+right+1 lines\n"
  "   -expanded        output medians and standard deviations instead of just\n"
  "                    averages\n"
  "   -header          useful particularly with the -expanded option to have a\n"
  "                    reminder of what column is what.\n"
  "   -cluster=k       cluster with k-means with given parameter (2-10 are best)\n"
  "   -cluster-sets=file.bed\n"
  "                    write out the original bed with the cluster label along\n"
  "                    with the accompanying region used for the clustering\n"
  "   -long-form=[label1,label2,...,labeln]\n"
  "                    output \"long form\" where each line is just the position\n"
  "                    and one of the values and the first column is the name of\n"
  "                    the file.\n"
  /* secret options:  they're not so important */
  /* "   -bed-ix=i        if a list of bed files is given, use only the ith one\n" */
  /* "   -bw-ix=i         if a list of bigWig files is given, only use the ith one\n" */
  );
}

struct agg_data
{
    int nrow;
    int ncol;
    int left;
    int right;
    int meta;
    int num_firsts;
    int num_seconds;
    int *indexes;
    char **first_names;
    char **second_names;
    double **data;
};

struct agg_data *init_agg_data(int left, int right, int meta, boolean firstbase, boolean nozero, int num_firsts,
			       int num_seconds, boolean expanded, struct slName *lf_labels)
/* init the output struct */
{
    struct agg_data *agg;
    struct slName *name;
    int i;
    int addone = 0;
    AllocVar(agg);
    agg->num_firsts = num_firsts;
    agg->num_seconds = num_seconds;
    agg->left = left;
    agg->right = right;
    agg->meta = meta;
    agg->nrow = left + right + meta;;
    agg->ncol = num_firsts * num_seconds;
    if (expanded)
	agg->ncol *= NUM_EXPANDED;
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
    AllocArray(agg->first_names, num_firsts);
    AllocArray(agg->second_names, num_seconds);
    i = 0;
    for (name = lf_labels; name != NULL; name = name->next)
    {
	if (i < num_firsts)
	    agg->first_names[i] = cloneString(name->name);
	else
	    agg->second_names[i-num_firsts] = cloneString(name->name);
	i++;
    }
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
    freeMem(agg->first_names);
    freeMem(agg->second_names);
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
		sum += pow(one_pbm_col[j] - mean,2);
	    if (size > 1)
		sd = sqrt(sum/(size-1));
	    else
		sd = na;
	    agg->data[i][offset] = mean;
	    if (expanded)
	    {
		agg->data[i][offset+1] = doubleMedian(size, one_pbm_col);
		agg->data[i][offset+2] = sd;
		agg->data[i][offset+3] = (double)size;
		agg->data[i][offset+4] = sum;
	    }
	}
	else
	{
	    agg->data[i][offset] = na;
	    if (expanded)
	    {
		agg->data[i][offset+1] = na;
		agg->data[i][offset+2] = na;
 		agg->data[i][offset+3] = 0;
 		agg->data[i][offset+4] = na;
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

static void output_data_of_line(FILE *out, struct agg_data *agg, boolean expanded, int i, int j)
{
    char buf[LONG_NUMBER];
    /* mean */
    NAorNumToS(buf, agg->data[i][j], -1);
    fprintf(out, "%s", buf);
    if (expanded)
    {
	double se = agg->data[i][j+2]/sqrt(agg->data[i][j+3]);
	double y_high = agg->data[i][j] + se;
	double y_low = agg->data[i][j] - se;
	/* median */
	NAorNumToS(buf, agg->data[i][j+1], -1);
	fprintf(out, "\t%s\t", buf);
	/* std dev */
	NAorNumToS(buf, agg->data[i][j+2], -1);
	fprintf(out, "%s\t", buf);
	/* number of data */
	fprintf(out, "%d\t", (int)agg->data[i][j+3]);
	/* sum squares */
	NAorNumToS(buf, agg->data[i][j+4], -1);
	fprintf(out, "%s\t", buf);
	/* SE */
	NAorNumToS(buf, se, -1);
	fprintf(out, "%s\t", buf);
	/* Y low */
	NAorNumToS(buf, y_low, -1);
	fprintf(out, "%s\t", buf);
	/* Y high */
	NAorNumToS(buf, y_high, -1);
	fprintf(out, "%s", buf);
    }
}

void output_agg_data(FILE *out, boolean expanded, boolean header, struct agg_data *agg, boolean long_form)
/* simply output the stuff */
{
    int i, j, k, l;
    char buf[LONG_NUMBER];
    if (long_form)
    /* currently there is no expanded form here */
    {
	if (header)
	{
	    fprintf(out, "Region\tSignal\tPosition\tMean");
	    if (expanded)
		fprintf(out, "\tMedian\tStd_Dev\tNum_Data\tSum_Squares\tStd_Err_Mean\tY_Low\tY_High\n");
	    else
		fprintf(out, "\n");
	}
	for (i = 0; i < agg->nrow; i++)
	{
	    k = 0;   /* index for first names (beds) */
	    l = 0;   /* index for second names (wigs) */
	    if ((agg->meta > 0) && (i == agg->left + agg->meta))
		fprintf(out, "# meta part ends here\n");
	    j = 0;
	    while (j < agg->ncol)
	    {
		fprintf(out, "%s\t%s\t%d\t", agg->first_names[k], agg->second_names[l], agg->indexes[i]);
		if (expanded)
		{
		    /* line out, newline  */
		    output_data_of_line(out, agg, expanded, i, j);
		    j += NUM_EXPANDED;
		}
		else
		{
		    output_data_of_line(out, agg, expanded, i, j);
		    j++;
		}
		fprintf(out, "\n");
		l++;
		if (l == agg->num_seconds)
		{
		    k++;
		    l = 0;
		}
	    }
	}
    }
    else
    {
	if (header)
	{
	    j = 0;
	    k = 1;
	    fprintf(out, "Position\t");
	    while (j < agg->ncol)
	    {
		if (expanded)
		{
		    fprintf(out, "Mean_%d\tMedian_%d\tStd_Dev_%d\tNum_Data_%d\tSum_Squares_%d\tStd_Err_Mean_%d\tY_High_%d\tY_Low_%d", k, k, k, k, k, k, k, k);
		    j += NUM_EXPANDED;
		}
		else
		{
		    fprintf(out, "Mean_%d", k);
		    j++;
		}
		fprintf(out, "%c", (j < agg->ncol) ? '\t' : '\n');
		k++;
	    }
	}
	for (i = 0; i < agg->nrow; i++)
	{
	    if ((agg->meta > 0) && (i == agg->left + agg->meta))
		fprintf(out, "# meta part ends here\n");
	    fprintf(out, "%d\t", agg->indexes[i]);
	    j = 0;
	    while (j < agg->ncol)
	    {
		/* line out, no newline */
		output_data_of_line(out, agg, expanded, i, j);
		j += (expanded) ? NUM_EXPANDED : 1;
		fprintf(out, "%c", (j < agg->ncol) ? '\t' : '\n');
	    }
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

static struct slName *setup_labels(char *long_form, boolean clustering, int k, struct slName *region_list, struct slName *wig_list,
			    struct slName **lf_labels_b, struct slName **lf_labels_w)
{
    int num_regions = slCount(region_list);
    int num_wigs = slCount(wig_list);
    struct slName *lf_labels = NULL;
    /* use the labels provided */
    if (long_form && !sameString(long_form, "on"))
    {
	lf_labels = slNameListFromComma(long_form);
	int num_labels = slCount(lf_labels);
	/* we should only have one or two labels if clustering is done */
	if (clustering)
	{
	    int i;
	    char buf[256];
	    /* handle the first case, where labels are only put on bed */
	    if (num_labels == 1)
		for (i = 1; i <= k; i++)
		{
		    safef(buf, sizeof(buf), "Cluster_%d", i);
		    slAddHead(&lf_labels, slNameNew(buf));
		}
	    /* the second case is to prefix the second label to each cluster name */
	    else if (num_labels == 2)
	    {
		char *second_name = lf_labels->next->name;
		struct slName *to_remove = slPopTail(&lf_labels);
		for (i = 1; i <= k; i++)
		{
		    safef(buf, sizeof(buf), "%s_cluster_%d", second_name, i);
		    slAddHead(&lf_labels, slNameNew(buf));
		}
		slNameFree(&to_remove);
	    }
	    else
		errAbort("Just use 0, 1, or 2 labels for -long-form when clustering");
	    slReverse(&lf_labels);
	}
	/* not clustering and labels are provided */
	else if (num_labels != num_regions + num_wigs)
	    errAbort("number of labels provided should equal the number of beds plus the number of bigWigs");
    }
    /* use labels from the filenames */
    else
    {
	if (*lf_labels_b || *lf_labels_w)
	{
	    if (*lf_labels_b && *lf_labels_w && (slCount(*lf_labels_b) == num_regions) && (slCount(*lf_labels_w) == num_wigs))
		lf_labels = slCat(*lf_labels_b, *lf_labels_w);
	    else
	    {
		if (*lf_labels_b)
		    slNameFreeList(lf_labels_b);
		if (*lf_labels_w)
		    slNameFreeList(lf_labels_w);
	    }
	}
	if (!lf_labels)
	{
	    struct slName *copy1 = slNameCloneList(region_list);
	    struct slName *copy2 = slNameCloneList(wig_list);
	    lf_labels = slCat(copy1, copy2);
	}
    }
    return lf_labels;
}

void bwtool_aggregate(struct hash *options, char *regions, unsigned decimals, double fill,
		      char *size_s, char *region_list_s, char *wig, char *tmp_dir, char *output_file)
/* aggregate - main */
{
    unsigned left = 0, right = 0;
    struct slName *region_list = slNameListFromComma(region_list_s);
    struct slName *wig_list = slNameListFromComma(wig);
    boolean firstbase = (hashFindVal(options, "firstbase") != NULL) ? TRUE : FALSE;
    boolean nozero = TRUE;
    int bed_ix = (int)sqlUnsigned((char *)hashOptionalVal(options, "bed-ix", "0"));
    int bw_ix = (int)sqlUnsigned((char *)hashOptionalVal(options, "bw-ix", "0"));
    boolean header = (hashFindVal(options, "header") != NULL) ? TRUE : FALSE;
    boolean use_start = (hashFindVal(options, "starts") != NULL) ? TRUE : FALSE;
    boolean use_end = (hashFindVal(options, "ends") != NULL) ? TRUE : FALSE;
    boolean expanded = (hashFindVal(options, "expanded") != NULL) ? TRUE : FALSE;
    boolean clustering = (hashFindVal(options, "cluster") != NULL) ? TRUE : FALSE;
    char *cluster_sets = (char *)hashFindVal(options, "cluster-sets");
    char *long_form = (char *)hashFindVal(options, "long-form");
    boolean do_long_form = (long_form != NULL);
    int meta = 0;
    struct slName *lf_labels = NULL, *lf_labels_b = NULL, *lf_labels_w = NULL;
    int k = (int)sqlUnsigned((char *)hashOptionalVal(options, "cluster", "0"));
    FILE *output;
    int num_regions = check_for_list_files(&region_list, &lf_labels_b, bed_ix);
    int num_wigs = check_for_list_files(&wig_list, &lf_labels_w, bw_ix);
    boolean mult_regions = (num_regions > 1);
    boolean mult_wigs = (num_wigs > 1);
    int num_parse = parse_left_right(size_s, &left, &right, &meta);
    boolean do_meta = (num_parse == 3);
    if (use_start && use_end)
	errAbort("cannot specify both -starts and -ends");
    if ((clustering) && ((k < 2) || (k > 10)))
	errAbort("k should be between 2 and 10\n");
    if ((mult_regions || mult_wigs) && clustering)
	errAbort("with clustering just specify one region list and one bigWig");
    if (do_meta && clustering)
	errAbort("at the moment clustering doesn't use meta\n");
    if (do_meta && (use_start || use_end || firstbase))
	errAbort("neither -firstbase nor -starts nor -ends are available with meta\n");
    if (firstbase)
	nozero = FALSE;
    if (mult_regions && mult_wigs)
	do_long_form = TRUE;
    lf_labels = setup_labels(long_form, clustering, k, region_list, wig_list, &lf_labels_b, &lf_labels_w);
    output = mustOpen(output_file, "w");
    if (!clustering)
    {
	int num_regions = slCount(region_list);
	int num_wigs = slCount(wig_list);
	struct slName *reg;
	struct slName *wig_name;
	struct agg_data *agg = NULL;
	struct metaBig *mbList = NULL;
	struct metaBig *mb;
	/* first calculate the meta if necessary as an average of all the regions in all the files */
	if (meta == -1)
	{
	    meta = calculate_meta_file_list(region_list);
	    fprintf(stderr, "calculated meta = %d bases\n", meta);
	}
	agg = init_agg_data(left, right, meta, firstbase, nozero, num_regions, num_wigs, expanded, lf_labels);
	for (wig_name = wig_list; wig_name != NULL; wig_name = wig_name->next)
	{
	    mb = metaBigOpenWithTmpDir(wig_name->name, tmp_dir, NULL);
	    slAddHead(&mbList, mb);
	}
	slReverse(&mbList);
	int offset = 0;
	if (!do_meta)
	{
	    for (reg = region_list; reg != NULL; reg = reg->next)
	    {
		struct bed6 *regions = load_and_recalculate_coords(reg->name, left, right, firstbase, use_start, use_end);
		struct slName *wig_name;
		for (mb = mbList; mb != NULL; mb = mb->next)
		{
		    struct perBaseMatrix *pbm = load_perBaseMatrix(mb, regions, fill);
		    do_summary(pbm, agg, expanded, offset);
		    offset += (expanded) ? NUM_EXPANDED : 1;
		    free_perBaseMatrix(&pbm);
		}
		bed6FreeList(&regions);
	    }
	}
	else
	/* the meta will be a fusion of three matrices */
	{
	    for (reg = region_list; reg != NULL; reg = reg->next)
	    {
		struct bed6 *regions_left = load_and_recalculate_coords(reg->name, left, 0, FALSE, TRUE, FALSE);
		struct bed6 *regions_right = load_and_recalculate_coords(reg->name, 0, right, FALSE, FALSE, TRUE);
		struct bed6 *regions_meta = (meta > 0) ? readBed6Soft(reg->name) : NULL;
		struct slName *wig_name;
		for (mb = mbList; mb != NULL; mb = mb->next)
		{
		    struct perBaseMatrix *pbm = load_perBaseMatrix(mb, regions_left, fill);
		    struct perBaseMatrix *right_pbm = load_perBaseMatrix(mb, regions_right, fill);
		    if (meta > 0)
		    {
			struct perBaseMatrix *meta_pbm = load_meta_perBaseMatrix(mb, regions_meta, meta, fill);
			fuse_pbm(&pbm, &meta_pbm, TRUE);
		    }
		    fuse_pbm(&pbm, &right_pbm, TRUE);
		    do_summary(pbm, agg, expanded, offset);
		    offset += (expanded) ? NUM_EXPANDED : 1;
		    free_perBaseMatrix(&pbm);
		}
		bed6FreeList(&regions_left);
		bed6FreeList(&regions_right);
		if (meta > 0)
		    bed6FreeList(&regions_meta);
	    }
	}
	while ((mb = slPopHead(&mbList)) != NULL)
	    metaBigClose(&mb);
    	output_agg_data(output, expanded, header, agg, do_long_form);
	free_agg_data(&agg);
    }
    else
    {
	struct agg_data *agg = init_agg_data(left, right, 0, firstbase, nozero, 1, k, FALSE, lf_labels);
	struct metaBig *mb = metaBigOpenWithTmpDir(wig_list->name, tmp_dir, NULL);
	struct bed6 *regions = load_and_recalculate_coords(region_list->name, left, right, firstbase, use_start, use_end);
	struct bed6 *orig_regions = readBed6Soft(region_list->name);
	struct perBaseMatrix *pbm = load_perBaseMatrix(mb, regions, fill);
	if (cluster_sets)
	    perBaseMatrixAddOrigRegions(pbm, orig_regions);
	struct cluster_bed_matrix *cbm = init_cbm_from_pbm(pbm, k);
	do_kmeans_sort(cbm, 0.001, TRUE);
	copy_centroids(cbm, agg);
	output_agg_data(output, FALSE, FALSE, agg, do_long_form);
	if (cluster_sets)
	    output_cluster_sets(cbm, cluster_sets);
	free_cbm(&cbm);
	metaBigClose(&mb);
	bed6FreeList(&regions);
	free_agg_data(&agg);
    }
    carefulClose(&output);
    slNameFreeList(&lf_labels);
    slNameFreeList(&region_list);
    slNameFreeList(&wig_list);
}
