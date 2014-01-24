#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/*****
 ** kmeans.c
 ** - a simple k-means clustering routine
 ** - returns the cluster labels of the data points in an array
 ** - here's an example
 **   extern int *k_means(double**, int, int, int, double, double**);
 **   ...
 **   int *c = k_means(data_points, num_points, dim, 20, 1e-4, 0);
 **   for (i = 0; i < num_points; i++) {
 **      printf("data point %d is in cluster %d\n", i, c[i]);
 **   }
 **   ...
 **   free(c);
 ** Parameters
 ** - array of data points (double **data)
 ** - number of data points (int n)
 ** - dimension (int m)
 ** - desired number of clusters (int k)
 ** - error tolerance (double t)
 **   - used as the stopping criterion, i.e. when the sum of
 **     squared euclidean distance (standard error for k-means)
 **     of an iteration is within the tolerable range from that
 **     of the previous iteration, the clusters are considered
 **     "stable", and the function returns
 **   - a suggested value would be 0.0001
 ** - output address for the final centroids (double **centroids)
 **   - user must make sure the memory is properly allocated, or
 **     pass the null pointer if not interested in the centroids
 ** References
 ** - J. MacQueen, "Some methods for classification and analysis
 **   of multivariate observations", Fifth Berkeley Symposium on
 **   Math Statistics and Probability, 281-297, 1967.
 ** - I.S. Dhillon and D.S. Modha, "A data-clustering algorithm
 **   on distributed memory multiprocessors",
 **   Large-Scale Parallel Data Mining, 245-260, 1999.
 ** Notes
 ** - this function is provided as is with no warranty.
 ** - the author is not responsible for any damage caused
 **   either directly or indirectly by using this function.
 ** - anybody is free to do whatever he/she wants with this
 **   function as long as this header section is preserved.
 ** Created on 2005-04-12 by
 ** - Roger Zhang (rogerz@cs.dal.ca)
 ** Modifications
 ** -
 ** Last compiled under Linux with gcc-3
 */

#include "common.h"
#include "obscure.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "basicBed.h"
#include "bigWig.h"
#include "metaBig.h"
#include "bigs.h"
#include "cluster.h"

static int perBaseWigJustLabelCmp(const void *a, const void *b)
/* for sorting after clustering */
{
    const struct perBaseWig *pbw_a = *((struct perBaseWig **)a);
    const struct perBaseWig *pbw_b = *((struct perBaseWig **)b);    
    int diff = pbw_a->label - pbw_b->label;
    return diff;
}

static int clear_na_rows(struct perBaseMatrix *pbm)
/* if an NA is encountered in the matrix row, set its label to -1 */
/* and put it at the beginning. return the number of rows like this */
{
    int i, j;
    int num_na = 0;
    for (i = 0; i < pbm->nrow; i++)
    {
	for (j = 0; j < pbm->ncol; j++)
	{
	    if (isnan(pbm->matrix[i][j]))
	    {
		pbm->array[i]->label = -1;
		pbm->array[i]->cent_distance = 0;
		num_na++;
		break;
	    }
	}
    }
    qsort(pbm->array, pbm->nrow, sizeof(pbm->array[0]), perBaseWigJustLabelCmp);
    for (i = 0; i < pbm->nrow; i++)
	pbm->matrix[i] = pbm->array[i]->data;
    return num_na;
}

struct cluster_bed_matrix *init_cbm_from_pbm(struct perBaseMatrix *pbm, int k)
/* initialize the cluster struct froma matrix */
{
    struct cluster_bed_matrix *cbm;
    int i;
    AllocVar(cbm);
    cbm->pbm = pbm;
    cbm->m = pbm->ncol;
    cbm->n = pbm->nrow;
    cbm->k = k;
    cbm->num_na = clear_na_rows(cbm->pbm);
    AllocArray(cbm->cluster_sizes, k);
    AllocArray(cbm->centroids, k);
    for (i = 0; i < k; i++)
    {
	AllocArray(cbm->centroids[i], cbm->m);
    }
    return cbm;
}

struct cluster_bed_matrix *init_cbm(struct metaBig *mb, struct bed6 *regions, int k, double fill)
/* initialize the cluster struct */
{
    struct perBaseMatrix *pbm = load_perBaseMatrix(mb, regions, fill);
    return init_cbm_from_pbm(pbm, k);
}

void free_cbm(struct cluster_bed_matrix **pCbm)
/* free up the cluster struct */
{
    struct cluster_bed_matrix *cbm = *pCbm;
    int i;
    free_perBaseMatrix(&cbm->pbm);
    freeMem(cbm->cluster_sizes);
    for (i = 0; i < cbm->k; i++)
	freeMem(cbm->centroids[i]);
    freeMem(cbm->centroids);
    freez(&cbm);
}

static int *k_means(struct cluster_bed_matrix *cbm, double t)
{
    /* output cluster label for each data point */
    int *labels; /* Labels for each cluster (size n) */
    int h, i, j; /* loop counters, of course :) */
    double old_error; 
    double error = DBL_MAX; /* sum of squared euclidean distance */
    double **tmp_centroids; /* centroids and temp centroids (size k x m) */
    int n = cbm->n;
    int m = cbm->m;
    int k = cbm->k;
    int pass = 1;
    AllocArray(labels, n);
    AllocArray(tmp_centroids, k);
    for (i = 0; i < k; i++)
	AllocArray(tmp_centroids[i], m);
    /* assert(data && k > 0 && k <= n && m > 0 && t >= 0); /\* for debugging *\/ */
    /* init ialization */
    for (i = 0, h = cbm->num_na; i < k; h += (cbm->n-cbm->num_na) / k, i++) 
    {
	/* pick k points as initial centroids */
	for (j = 0; j < m; j++)
	    cbm->centroids[i][j] = cbm->pbm->matrix[h][j];
    }
    /* main loop */
    do 
    {
	/* save error from last step */
	old_error = error;
	error = 0;
	/* clear old cbm->cluster_sizes and temp centroids */
	for (i = 0; i < k; i++)
	{
	    cbm->cluster_sizes[i] = 0;
	    for (j = 0; j < m; j++)
		tmp_centroids[i][j] = 0;
	}
	for (h = cbm->num_na; h < n; h++) 
	{
	    /* identify the closest cluster */
	    double min_distance = DBL_MAX;
	    for (i = 0; i < k; i++) 
	    {
		double distance = 0;
		for (j = 0; j < m; j++)
		    distance += pow(cbm->pbm->matrix[h][j] - cbm->centroids[i][j], 2);
		distance = sqrt(distance);
		if (distance < min_distance) 
		{
		    labels[h] = i;
		    min_distance = distance;
		}
	    }
	    /* update size and temp centroid of the destination cluster */
	    for (j = 0; j < m; j++)
		tmp_centroids[labels[h]][j] += cbm->pbm->matrix[h][j];
	    cbm->cluster_sizes[labels[h]]++;
	    /* update standard error */
	    cbm->pbm->array[h]->cent_distance = min_distance;
	    error += min_distance;
	}
	/* update all centroids */
	for (i = 0; i < k; i++) 
	    for (j = 0; j < m; j++) 
		cbm->centroids[i][j] = (cbm->cluster_sizes[i] > 0)  
		    ? (tmp_centroids[i][j] / cbm->cluster_sizes[i]) : 0;
    } while (fabs(error - old_error) > t);
    /* housekeeping */
    for (i = 0; i < k; i++)
	freeMem(tmp_centroids[i]);
    freeMem(tmp_centroids);
    return labels;
}

int sortthesize(const void *a, const void *b)
{
    int *row_a = *(int **)a;
    int *row_b = *(int **)b;
    int diff = row_b[0] - row_a[0];
    if (diff == 0)
	return row_a[1] - row_b[1];
    return diff;
}

int sortoriglabel(const void *a, const void *b)
{
    int *row_a = *(int **)a;
    int *row_b = *(int **)b;
    return row_a[1] - row_b[1];
}

static void exchange_labels(int *cluster_sizes, int k, int *labels, int start, int size)
/* rename labels based on size */
{
    int **array;
    int i;
    AllocArray(array, k);
    for (i = 0; i < k; i++)
    {
	AllocArray(array[i], 3);
	array[i][0] = cluster_sizes[i];
	array[i][1] = i;
    }
    qsort(array, k, sizeof(int *), sortthesize);
    for (i = 0; i < k; i++)
    {
	array[i][2] = i;
	cluster_sizes[i] = array[i][0];
    }
    qsort(array, k, sizeof(int *), sortoriglabel);
    for (i = start; i < size; i++)
	labels[i] = array[labels[i]][2];
    for (i = 0; i < k; i++)
	freeMem(array[i]);
    freeMem(array);
}

void do_kmeans_sort(struct cluster_bed_matrix *cbm, double t, boolean sort)
/* clusters but also sorts the labels by cluster size */
{
    int i = 0;
    int *labels = k_means(cbm, t);
    /* if (sort) */
    /* 	exchange_labels(cbm->cluster_sizes, cbm->k, labels, cbm->num_na, cbm->pbm->nrow); */
    for (i = cbm->num_na; i < cbm->pbm->nrow; i++)
	cbm->pbm->array[i]->label = labels[i];
    qsort(cbm->pbm->array, cbm->pbm->nrow, sizeof(cbm->pbm->array[0]), perBaseWigLabelCmp);
    for (i = 0; i < cbm->pbm->nrow; i++)
	cbm->pbm->matrix[i] = cbm->pbm->array[i]->data;
    freeMem(labels);
}

void do_kmeans(struct cluster_bed_matrix *cbm, double t)
/* the main clustering function.  labels matrix rows and reorders */
{
    return do_kmeans_sort(cbm, t, FALSE);
}
