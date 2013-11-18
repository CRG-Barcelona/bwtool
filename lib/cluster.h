#ifndef CLUSTER_H
#define CLUSTER_H

struct cluster_bed_matrix
/* struct holding all the needed info for k-means clustering */
{
    struct perBaseMatrix *pbm;
    struct bed *orig_list;
    int n;
    int m;
    int k;
    int num_na;
    int *cluster_sizes;
    double **centroids;
};

struct cluster_bed_matrix *init_cbm_from_pbm(struct perBaseMatrix *pbm, int k);
/* initialize the cluster struct froma matrix */

struct cluster_bed_matrix *init_cbm(struct metaBig *mb, struct bed6 *regions, int k, double fill);
/* initialize the cluster struct */

void free_cbm(struct cluster_bed_matrix **pCbm);
/* free up the cluster struct */

void do_kmeans_sort(struct cluster_bed_matrix *cbm, double t, boolean sort);
/* clusters but also sorts the labels by cluster size */

void do_kmeans(struct cluster_bed_matrix *cbm, double t);
/* the main clustering function.  labels matrix rows and reorders */

#endif /* cluster.h */
