/* bwtool_matrix - with a bed file of regions, select  */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdint.h>
#include <jkweb/common.h>
#include <jkweb/obscure.h>
#include <jkweb/linefile.h>
#include <jkweb/hash.h>
#include <jkweb/options.h>
#include <jkweb/sqlNum.h>
#include <jkweb/basicBed.h>
#include <jkweb/bigWig.h>
#include <beato/bigs.h>
#include <beato/metaBig.h>
#include "bwtool.h"
#include <beato/cluster.h>
#include "bwtool_shared.h"

/* bwtool binary matrix magic bytes:
     BwTool\x91\x90
     42 77 54 6f 6f 6c 91 90
*/
const uint64_t BWTOOL_MAGIC = 0x90916c6f6f547742;

/* bwtool binary matrix file format version number. */
const uint64_t BWTOOL_FILEFORMAT_VERSION = 1;

/* bwtool binary matrix type. */
const uint64_t BWTOOL_BINARY_MATRIX = 0;
/* bwtool binary cluster matrix type. */
const uint64_t BWTOOL_BINARY_CLUSTER_MATRIX = 1;

struct binary_matrix_header
{
    /* Magic bytes: BWTOOL_MAGIC */
    uint64_t bwtool_magic;
    /* File format version: BWTOOL_FILEFORMAT_VERSION */
    uint64_t bwtool_fileformat_version;
    /* bwtool binary (cluster) matrix type: BWTOOL_BINARY_MATRIX or BWTOOL_BINARY_CLUSTER_MATRIX */
    uint64_t bwtool_matrix_type;
    /* Number of rows in the matrix. */
    uint64_t nbr_rows;
    /* Number of columns in the matrix. */
    uint64_t nbr_columns;
    /* Size of float (double) in bytes. */
    uint64_t float_byte_size;
    /* Matrix length in bytes: nbr_rows * nbr_columns * sizeof(double). */
    uint64_t matrix_length_in_bytes;
    /* Offset of binary matrix: 4096 */
    uint64_t matrix_offset;
    /* Offset of text metadata file: multiple of 4096 (if 0: no text metadata file). */
    uint64_t text_metadata_offset;
};

struct binary_matrix_header make_binary_matrix_header(struct perBaseMatrix *pbm, boolean keep_bed)
{
    uint64_t nbr_rows = (uint64_t) pbm->nrow;
    uint64_t nbr_columns = (uint64_t) pbm->ncol;
    uint64_t float_byte_size = (uint64_t) sizeof(double);
    uint64_t matrix_length_in_bytes = nbr_rows * nbr_columns * sizeof(double);
    uint64_t matrix_offset = 4096;
    uint64_t text_metadata_offset = 0;

    if (keep_bed) {
	if ( ((matrix_offset + matrix_length_in_bytes) % 4096) != 0)
	{
	    text_metadata_offset = ( ( (matrix_offset + matrix_length_in_bytes) / 4096 ) + 1 ) * 4096;
	}
	else
	{
	    text_metadata_offset = matrix_offset + matrix_length_in_bytes;
	}
    }

    struct binary_matrix_header binary_matrix_header_instance =
    {
	BWTOOL_MAGIC,
	BWTOOL_FILEFORMAT_VERSION,
	BWTOOL_BINARY_MATRIX,
	nbr_rows,
	nbr_columns,
	float_byte_size,
	matrix_length_in_bytes,
	matrix_offset,
	text_metadata_offset,
    };

    return binary_matrix_header_instance;
}

struct binary_matrix_header make_binary_cluster_matrix_header(struct cluster_bed_matrix *cbm, boolean keep_bed)
{
    uint64_t nbr_rows = (uint64_t) cbm->pbm->nrow;
    uint64_t nbr_columns = (uint64_t) cbm->pbm->ncol;
    uint64_t float_byte_size = (uint64_t) sizeof(double);
    uint64_t matrix_length_in_bytes = nbr_rows * nbr_columns * sizeof(double);
    uint64_t matrix_offset = 4096;
    uint64_t text_metadata_offset = 0;

    if (keep_bed) {
	if ( ((matrix_offset + matrix_length_in_bytes) % 4096) != 0)
	{
	    text_metadata_offset = ( ( (matrix_offset + matrix_length_in_bytes) / 4096 ) + 1 ) * 4096;
	}
	else
	{
	    text_metadata_offset = matrix_offset + matrix_length_in_bytes;
	}
    }

    struct binary_matrix_header binary_matrix_header_instance =
    {
	BWTOOL_MAGIC,
	BWTOOL_FILEFORMAT_VERSION,
	BWTOOL_BINARY_CLUSTER_MATRIX,
	nbr_rows,
	nbr_columns,
	float_byte_size,
	matrix_length_in_bytes,
	matrix_offset,
	text_metadata_offset,
    };

    return binary_matrix_header_instance;
}

void usage_matrix()
/* Explain usage of matrix program and exit. */
{
errAbort(
  "bwtool matrix - extract data and output in a tab-delimited way that is\n"
  "   easily used as a matrix by other programs.\n"
  "usage:\n"
  "   bwtool matrix up:down regions.bed input1.bw,input2... output.txt\n"
  "   bwtool matrix up:down regions.bed bigWigs.lst output.txt\n\n"
  "   bwtool matrix up:meta:down regions.bed input1.bw,input2... output.txt\n"
  "   bwtool matrix up:meta:down regions.bed bigWigs.lst output.txt\n\n"
  "If more than one bigWig is specified, then the matrix from the second\n"
  "matrix is fused to the first, and the third to the second, etc. in the\n"
  "same left-to-right order as the comma-list of bigWigs\n\n"
  "options:\n"
  "   -binary-matrix  write a binary matrix file instead of a tab-delimited\n"
  "                   file which can be imported faster by programs that support\n"
  "                   reading of binary files\n"
  "   -keep-bed       in this case output the original bed loci in the first\n"
  "                   columns of the output and output data as comma-separated\n"
  "   -starts         use starts of bed regions as opposed to the middles\n"
  "   -ends           use ends of bed regions as opposed to the middles\n"
  "   -tiled-averages=n\n"
  "                   break the up:down sized region into regions of n bases\n"
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

void na_or_num(FILE *out, double num, int decimals)
/* repeated function testing NA-ness of the number then outputting. */
{
    if (isnan(num))
	fprintf(out, "NA");
    else
	fprintf(out, "%0.*f", decimals, num);
}

void output_centroids(struct cluster_bed_matrix *cbm, char *centroid_file, int decimals)
/* strictly the centroids */
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
/* non-long output */
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
	{
	    na_or_num(out, pbw->data[j], decimals);
	    fprintf(out, "%c", (j == pbw->len-1) ? '\n' : '\t');
	}
    }
    carefulClose(&out);
}

void output_binary_cluster_matrix(struct cluster_bed_matrix *cbm, boolean keep_bed, char *outputfile)
/* Binary cluster matrix output */
{
    FILE *out = mustOpen(outputfile, "w");
    int i;

    /* Write the header. */
    struct binary_matrix_header binary_matrix_header_instance = make_binary_cluster_matrix_header(cbm, keep_bed);
    fwrite(&binary_matrix_header_instance, sizeof(binary_matrix_header_instance), 1, out);

    /* Put binary cluster matrix at offset 4096. */
    fseek(out, 4096, SEEK_SET);

    /* Write the binary cluster matrix. */
    for (i = 0; i < cbm->pbm->nrow; i++)
    {
	struct perBaseWig *pbw = cbm->pbm->array[i];
	fwrite(pbw->data, sizeof(double), pbw->len, out);
    }

    /* Put the text metadata at a multiple of 4096. */
    fseek(out, binary_matrix_header_instance.text_metadata_offset, SEEK_SET);

    /* Write the text metadata after the matrix. */
    for (i = 0; i < cbm->pbm->nrow; i++)
    {
	struct perBaseWig *pbw = cbm->pbm->array[i];
	if (keep_bed)
	{
	    fprintf(out, "%s\t%d\t%d\t%s\t%d\t%c\t", pbw->chrom, pbw->chromStart, pbw->chromEnd, pbw->name, pbw->score, pbw->strand[0]);
	}
	fprintf(out, "%d\t", pbw->label);
	fprintf(out, "%f\n", pbw->cent_distance);
    }
    carefulClose(&out);
}

void output_cluster_matrix_long(struct cluster_bed_matrix *cbm, struct slName *labels, boolean keep_bed, char *outputfile, boolean header)
/* For handling long-form cluster output */
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

void output_matrix_long(struct perBaseMatrix *pbm, int decimals, struct slName *labels, boolean keep_bed, int left,
			int right, int tile, boolean header, char *outputfile)
/* long output.  right this is just patching things up.  this and some other stuff could be combined */
/* with aggregate some day. */
{
    FILE *out = mustOpen(outputfile, "w");
    int i,j,k, lr_pos;
    int n_labels = slCount(labels);
    int unfused_cols = pbm->ncol / n_labels;
    assert(unfused_cols == (left+right)/tile);
    struct slName *lab;
    if (header)
    {
	if (keep_bed)
	    fprintf(out, "chrom\tchromStart\tchromEnd\tname\tscore\tstrand\tSignal\tPosition\tValue\n");
	else
	    fprintf(out, "Signal\tRegion\tPosition\tValue\n");
    }
    /* Do label, region, Position, Value */
    for (lr_pos = -1 *left, k = 0; (lr_pos <= right) && (k < unfused_cols); lr_pos += (lr_pos + tile == 0) ? 2*tile : tile, k++)
    {
	for (lab = labels, j = 0; (lab != NULL) && (j < n_labels); lab = lab->next, j++)
	{
	    for (i = 0; i < pbm->nrow; i++)
	    {
		int jj = j*unfused_cols + k;
		struct perBaseWig *pbw = pbm->array[i];
		if (keep_bed)
		{
		    char strand = pbw->strand[0];
		    char *chrom = pbw->chrom;
		    int chromStart = k*tile + pbw->chromStart;
		    if (strand == '-')
			chromStart = pbw->chromEnd - k*tile - 1;
		    fprintf(out, "%s\t%d\t%d\t", chrom, chromStart, chromStart + tile);
		    fprintf(out, "%s\t", pbw->name);
		    fprintf(out, "%d\t%c\t", pbw->score, strand);
		}
		fprintf(out, "%s\t", lab->name);
		if (!keep_bed)
		    fprintf(out, "%s\t", (pbw->name) ? pbw->name : ".");
		fprintf(out, "%d\t", lr_pos);
		na_or_num(out, pbm->matrix[i][jj], decimals);
		fprintf(out, "\n");
	    }
	}
    }
    carefulClose(&out);
}

void output_matrix(struct perBaseMatrix *pbm, int decimals, boolean keep_bed, char *outputfile)
/* the simplest output */
{
    FILE *out = mustOpen(outputfile, "w");
    int i,j;
    for (i = 0; i < pbm->nrow; i++)
    {
	struct perBaseWig *pbw = pbm->array[i];
	if (keep_bed)
	    fprintf(out, "%s\t%d\t%d\t%s\t%d\t%c\t", pbw->chrom, pbw->chromStart, pbw->chromEnd, pbw->name, pbw->score, pbw->strand[0]);
	for (j = 0; j < pbw->len; j++)
	{
	    na_or_num(out, pbw->data[j], decimals);
	    fprintf(out, "%c", (j == pbw->len-1) ? '\n' : '\t');
	}
    }
    carefulClose(&out);
}

void output_binary_matrix(struct perBaseMatrix *pbm, boolean keep_bed, char *outputfile)
/* Binary matrix output */
{
    FILE *out = mustOpen(outputfile, "w");
    int i;

    /* Write the header. */
    struct binary_matrix_header binary_matrix_header_instance = make_binary_matrix_header(pbm, keep_bed);
    fwrite(&binary_matrix_header_instance, sizeof(binary_matrix_header_instance), 1, out);

    /* Put matrix at offset 4096. */
    fseek(out, 4096, SEEK_SET);

    /* Write the binary matrix. */
    for (i = 0; i < pbm->nrow; i++)
    {
	struct perBaseWig *pbw = pbm->array[i];
        fwrite(pbw->data, sizeof(double), pbw->len, out);
    }

    /* Write the BED file after the matrix, if requested. */
    if (keep_bed){
        /* Put BED file at a multiple of 4096. */
        fseek(out, binary_matrix_header_instance.text_metadata_offset, SEEK_SET);

	for (i = 0; i < pbm->nrow; i++)
	{
	    struct perBaseWig *pbw = pbm->array[i];
	    fprintf(out, "%s\t%d\t%d\t%s\t%d\t%c\n", pbw->chrom, pbw->chromStart, pbw->chromEnd, pbw->name, pbw->score, pbw->strand[0]);
	}
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

void check_list(struct bed6 *list)
{
    struct bed6 *one;
    for (one = list; one != NULL; one = one->next)
	uglyf("   %s\t%d\t%d\n", one->chrom, one->chromStart, one->chromEnd);
}

void bwtool_matrix(struct hash *options, char *favorites, char *regions, unsigned decimals,
		   double fill, char *range_s, char *bigfile, char *tmp_dir, char *outputfile)
/* bwtool_matrix - main for matrix-creation program */
{
    boolean do_k = (hashFindVal(options, "cluster") != NULL) ? TRUE : FALSE;
    boolean do_tile = (hashFindVal(options, "tiled-averages") != NULL) ? TRUE : FALSE;
    boolean do_binary_matrix = (hashFindVal(options, "binary-matrix") != NULL) ? TRUE : FALSE;
    boolean keep_bed = (hashFindVal(options, "keep-bed") != NULL) ? TRUE : FALSE;
    boolean starts = (hashFindVal(options, "starts") != NULL) ? TRUE : FALSE;
    boolean ends = (hashFindVal(options, "ends") != NULL) ? TRUE : FALSE;
    boolean lf_header = (hashFindVal(options, "long-form-header") != NULL) ? TRUE : FALSE;
    char *centroid_file = (char *)hashFindVal(options, "cluster-centroids");
    char *long_form = (char *)hashFindVal(options, "long-form");
    boolean do_long_form = (long_form != NULL);
    unsigned left = 0, right = 0;
    int meta = 0;
    int num_parse = parse_left_right(range_s, &left, &right, &meta);
    boolean do_meta = (num_parse == 3);
    int k = (int)sqlUnsigned((char *)hashOptionalVal(options, "cluster", "0"));
    int tile = (int)sqlUnsigned((char *)hashOptionalVal(options, "tiled-averages", "1"));
    if ((do_k) && ((k < 2) || (k > 10)))
	errAbort("k should be between 2 and 10\n");
    if ((do_tile) && (tile < 2))
	errAbort("tiling should be done for larger regions");
    if ((left % tile != 0) || (right % tile != 0))
	errAbort("tiling should be multiple of both up and down values");
    if (do_meta && starts && ends)
	warn("meta uses -starts and -ends anyway");
    else if ((do_meta) && (starts || ends))
	warn("-starts and -ends both automatically used with meta");
    if (do_meta && do_tile)
	errAbort("meta not compatible with -tile... yet");
    if (do_binary_matrix && do_long_form)
	errAbort("Writing binary matrix is not compatible with -long-form... yet");
    struct slName *bw_names = slNameListFromComma(bigfile);
    struct slName *bw_name;
    struct slName *labels_from_file = NULL;
    int num_bigwigs = check_for_list_files(&bw_names, &labels_from_file, 0);
    struct slName *labels = setup_labels(long_form, bw_names, &labels_from_file);
    struct bed6 *regs = NULL;
    struct bed6 *regions_left = NULL, *regions_right = NULL, *regions_meta = NULL;
    struct perBaseMatrix *pbm = NULL;
    int i;
    if (do_meta)
    {
	if (meta == -1)
	{
	    meta = calculate_meta_file(regions);
	    fprintf(stderr, "calculated meta = %d bases\n", meta);
	}
	regions_left = load_and_recalculate_coords(regions, left, 0, FALSE, TRUE, FALSE);
	regions_right = load_and_recalculate_coords(regions, 0, right, FALSE, FALSE, TRUE);
	regions_meta = (meta > 0) ? readBed6Soft(regions) : NULL;
    }
    else
	regs = load_and_recalculate_coords(regions, left, right, FALSE, starts, ends);
    for (bw_name = bw_names; bw_name != NULL; bw_name = bw_name->next)
    {
	struct metaBig *mb = metaBigOpenWithTmpDir(bw_name->name, tmp_dir, NULL);
	if (do_meta)
	{
	    struct perBaseMatrix *one_pbm = load_perBaseMatrix(mb, regions_left, fill);
	    struct perBaseMatrix *right_pbm = load_perBaseMatrix(mb, regions_right, fill);
	    if (meta > 0)
	    {
		struct perBaseMatrix *meta_pbm = load_meta_perBaseMatrix(mb, regions_meta, meta, fill);
		fuse_pbm(&one_pbm, &meta_pbm, TRUE);
	    }
	    fuse_pbm(&one_pbm, &right_pbm, TRUE);
	    fuse_pbm(&pbm, &one_pbm, FALSE);
	}
	else
	{
	    struct perBaseMatrix *one_pbm = (do_tile) ? load_ave_perBaseMatrix(mb, regs, tile, fill) :
		load_perBaseMatrix(mb, regs, fill);
	    fuse_pbm(&pbm, &one_pbm, FALSE);
	}
	metaBigClose(&mb);
    }
    if (do_k)
    {
	struct cluster_bed_matrix *cbm = NULL;
	/* ordered by cluster with label in first column */
	cbm = init_cbm_from_pbm(pbm, k);
	do_kmeans_sort(cbm, 0.001, TRUE);
	if (do_long_form)
	{
	    output_cluster_matrix_long(cbm, labels, keep_bed, outputfile, lf_header);
	}
	else
	{
	    if (do_binary_matrix)
	    {
		output_binary_cluster_matrix(cbm, keep_bed, outputfile);
	    }
	    else
	    {
		output_cluster_matrix(cbm, decimals, keep_bed, outputfile);
	    }
	}
	if (centroid_file)
	{
	    output_centroids(cbm, centroid_file, decimals);
	}
	free_cbm(&cbm);
    }
    else
    {
	if (do_long_form)
	{
	    output_matrix_long(pbm, decimals, labels, keep_bed, left, right, tile, lf_header, outputfile);
	}
	else
	{
	    if (do_binary_matrix)
	    {
		output_binary_matrix(pbm, keep_bed, outputfile);
	    }
	    else
	    {
		output_matrix(pbm, decimals, keep_bed, outputfile);
	    }
	}
	/* unordered, no label  */
	free_perBaseMatrix(&pbm);
    }
    if (do_meta)
    {
	bed6FreeList(&regions_left);
	bed6FreeList(&regions_right);
	if (meta > 0)
	    bed6FreeList(&regions_meta);
    }
    else
	bed6FreeList(&regs);
}
