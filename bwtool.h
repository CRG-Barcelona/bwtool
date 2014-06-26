#ifndef BWTOOL_H
#define BWTOOL_H

#ifndef FLOAT_H
#include <float.h>
#endif

enum bw_op_type
{
    invalid = 0,
    not_equal = 1,
    equal = 2,
    less = 3,
    less_equal = 4,
    more = 5,
    more_equal = 6,
    mask = 7,
};

enum bw_op_type get_bw_op_type(char *thresh_type, boolean inverse);
/* speed things up inside loop.  we don't want a string comparison for every datapoint */
/* to know what kind of operation to perform. */

void usage_remove();
/* Explain removal usage and exit. */

void usage_fill();
/* Explain fill usage and exit. */

void usage_shift();
/* Explain usage and exit. */

void usage_find();
/* Explain find usage and exit. */

void usage_distrib();
/* Explain usage of distribution program and exit. */

void usage_chromgraph();
/* Explain usage of chromgraph program and exit. */

void usage_distance();
/* Explain usage of distance calculation program and exit. */

void usage_aggregate();
/* Explain usage of aggregate program and exit. */

void usage_random();
/* Explain usage of randomized data program and exit. */

void usage_summary();
/* Explain usage of the summarize program and exit. */

void usage_matrix();
/* Explain usage of the matrix-creation program and exit. */

void usage_roll();
/* Explain the usage of the rolling-mean program and exit. */

void usage_sax();
/* Explain usage of the sax program and exit. */

void usage_split();
/* Explain usage of the splitting program and exit. */

void usage_lift();
/* Explain usage of the lift program */

void usage_paste();
/* Explain usage of paste program and exit. */

void usage_window();
/* Explain usage of the window-tiling program and exit */

void usage_extract();
/* Explain usage and exit. */

void bwtool_remove(struct hash *options, char *favorites, char *regions, unsigned decimals, enum wigOutType wot,
		   boolean condense, boolean wig_only, char *thresh_type, char *val_or_file, char *bigfile, char *tmp_dir,
		   char *outputfile);
/* bwtool_remove - main for removal program */

void bwtool_fill(struct hash *options, char *favorites, char *regions, unsigned decimals, enum wigOutType wot,
		 boolean condense, char *val_s, char *bigfile, char *tmp_dir, char *outputfile);
/* bwtool_fill - main for filling program */

void bwtool_find_extrema(struct hash *options, char *favorites, char *regions, unsigned decimals,
			 double fill, char *bigfile, char *tmp_dir, char *outputfile);
/* find local extrema */

void bwtool_find_thresh(struct hash *options, char *favorites, char *regions, double fill,
			char *thresh_type, char *thresh_s, char *bigfile, char *tmp_dir, char *outputfile);
/* find regions fitting a specified threshold */

void bwtool_find_max(struct hash *options, char *favorites, char *regions, double fill,
		     char *bigfile, char *tmp_dir, char *outputfile);
/* find max points in a range */

void bwtool_distrib(struct hash *options, char *favorites, char *regions, unsigned decimals,
		    char *bigfile, char *tmp_dir, char *outputfile);
/* make plotting information about the bigWig's distribution */

void bwtool_chromgraph(struct hash *options, char *favorites, char *regions, unsigned decimals,
		       double fill, char *bigfile, char *tmp_dir, char *outputfile);
/* bwtool_chromgraph - main for making the chromgraph file */

void bwtool_aggregate(struct hash *options, char *regions, unsigned decimals, double fill,
		      char *size_s, char *region_list_s, char *wig, char *tmp_dir, char *output_file);
/* aggregate - main */

void bwtool_random(struct hash *options, char *favorites, char *regions, unsigned decimals,
		   double fill, char *num_s, char *size_s, char *bigfile, char *tmp_dir, char *output_file);
/* random - main */

void bwtool_summary(struct hash *options, char *favorites, char *regions, unsigned decimals,
		    double fill, char *loci_s, char *bigfile, char *tmp_dir, char *outputfile);
/* bwtool_summary - main for the summarize program */

void bwtool_shift(struct hash *options, char *favorites, char *regions, unsigned decimals, enum wigOutType wot,
		  boolean condense, char *val_s, char *bigfile, char *tmp_dir, char *outputfile);
/* bwtool_shift - main for shifting program */

void bwtool_split(struct hash *options, char *regions, char *size_s, char *bigfile, char *tmp_dir, char *outputfile);
/* bwtool_split - main for the splitting program */

void bwtool_matrix(struct hash *options, char *favorites, char *regions, unsigned decimals,
		   double fill, char *range_s, char *bigfile, char *tmp_dir, char *outputfile);
/* bwtool_matrix - main for matrix-creation program */

void bwtool_autocorr(struct hash *options, char *favorites, char *regions, unsigned decimals,
		     double fill, char *bigfile, char *tmp_dir, char *outputfile);
/* bwtool_autocorr - main for autocorrelation program */

void bwtool_sax(struct hash *options, char *favorites, char *regions, unsigned decimals,
		char *alpha_s, char *bigfile, char *tmp_dir, char *outputfile);
/* bwtool_sax - main for the sax symbol program */

void bwtool_lift(struct hash *options, char *favorites, char *regions, unsigned decimals,
		 enum wigOutType wot, char *bigfile, char *tmp_dir, char *chainfile, char *outputfile);
/* bwtool_lift - main for lifting program */

void bwtool_paste(struct hash *options, char *favorites, char *regions, unsigned decimals, double fill,
		  enum wigOutType wot, struct slName **p_files, char *tmp_dir, char *output_file);
/* bwtool_paste - main for paste program */

void bwtool_window(struct hash *options, char *favorites, char *regions, unsigned decimals, double fill,
                   char *size_s, char *bigfile, char *tmp_dir, char *outputfile);
/* bwtool_window - main for the windowing program */

void bwtool_extract(struct hash *options, char *regions, unsigned decimals, double fill,
		    char *style, char *bigfile, char *tmp_dir, char *output_file);
/* bwtool_extract - main for the extract program */

void bwtool_roll(struct hash *options, char *favorites, char *regions, unsigned decimals, double fill,
		 enum wigOutType wot, char *command, char *size_s, char *bigfile, char *tmp_dir, char *outputfile);
/* bwtool_roll - main for the rolling-mean program */

#endif /* BWTOOL_H */
