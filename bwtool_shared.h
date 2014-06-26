#ifndef BWTOOL_SHARED_H
#define BWTOOL_SHARED_H

#include <jkweb/common.h>
#include <beato/bigs.h>

struct bed6 *load_and_recalculate_coords(char *list_file, int left, int right, boolean firstbase, boolean starts, boolean ends);
/* do the coordinate recalculation */

int check_for_list_files(struct slName **pList, struct slName **lf_list_labels, int ix);
/* check if the comma-lists contain actual files or if they contain files listing files */
/* either way, it should expand the file into a list and return the size of the list */
/* if it's the case that the comma-list are files, then nothing should change. */

int parse_left_right(char *size_s, unsigned *pleft, unsigned *pright, int *pmeta);
/* parse the "left:right" or "left:meta:right" from the command */
/* return the number of args: 2 or 3 */

void writeBw(char *inName, char *outName, struct hash *chromSizeHash);
/* shared func */

struct metaBig *metaBigOpen_check(char *bigfile, char *tmp_dir, char *regions);
/* A wrapper for metaBigOpen that does some checking and erroring */

void fuse_pbm(struct perBaseMatrix **pBig, struct perBaseMatrix **pTo_add, boolean add_coords);
/* not this makes perhaps-illegal perBaseWigs where the chromEnd-chromStart are not the */
/* same as the len... which may break things somewhere if this were ever library-ized */

int calculate_meta_file(char *file_name);
/* from all the beds one region file, get a single average */

int calculate_meta_file_list(struct slName *region_list);
/* from all the beds in all the region files, get a single average */

#endif /* BWTOOL_SHARED_H */
