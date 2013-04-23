#ifndef FAVS_H
#define FAVS_H

#ifndef HASH_H
#include "hash.h"
#endif

#define FAVS_FILE "~/.favorites.txt"

struct hash *favs_load_hash(char *filename);
/* load up the favorites hash from the file. NULL filename loads the default. */

#endif
