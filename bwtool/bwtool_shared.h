#ifndef BWTOOL_SHARED_H
#define BWTOOL_SHARED_H

#include "common.h"
#include "bigs.h"

struct bed6 *load_and_recalculate_coords(char *list_file, int left, int right, boolean firstbase, boolean starts, boolean ends);
/* do the coordinate recalculation */

int check_for_list_files(struct slName **pList, struct slName **lf_list_labels);
/* check if the comma-lists contain actual files or if they contain files listing files */
/* either way, it should expand the file into a list and return the size of the list */
/* if it's the case that the comma-list are files, then nothing should change. */

void parse_left_right(char *size_s, unsigned *pleft, unsigned *pright);
/* parse the "left:right" from the command */

void writeBw(char *inName, char *outName, struct hash *chromSizeHash);
/* shared func */

#endif /* BWTOOL_SHARED_H */
