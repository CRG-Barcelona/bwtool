#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* library that provides some functionality to the utilities for keeping track of */
/* files and URLs easily */

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "favs.h"

struct hash *favs_load_hash(char *filename)
/* load up the favorites hash from the file. NULL filename loads the default. */
{
    struct lineFile *lf = NULL;
    char file[1024];
    struct hash *hash = NULL;
    char *words[2];
    if (filename == NULL)
    {
	char *home_var = getenv("HOME");
	if (home_var) 
	    safef(file, sizeof(file), "%s/.favorites.txt", home_var);
	else
	    safef(file, sizeof(file), "favorites.txt");;
    }
    else
	safef(file, sizeof(file), "%s", filename);
    lf = lineFileOpen(file, TRUE);
    hash = newHash(20);
    while (lineFileRowTab(lf, words))
	hashAdd(hash, words[0], words[1]);
    return hash;
}
