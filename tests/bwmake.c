/* bwmake - doesn't do anything except make a bigWig */
/*   ** this isn't a full-featured program.  I'ts meant */
/*   ** just for the test script. */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif 

#include <jkweb/common.h>
#include <jkweb/hash.h>
#include <jkweb/linefile.h>
#include <jkweb/sqlNum.h>
#include <jkweb/bigWig.h>
#include <jkweb/bwgInternal.h>

void writeBw(char *inName, char *outName, struct hash *chromSizeHash)
/* copied from bwtool_shared.c */
{
    struct lm *lm = lmInit(0);
    struct bwgSection *sectionList = bwgParseWig(inName, TRUE, chromSizeHash, 1024, lm);
    if (sectionList == NULL)
	errAbort("%s is empty of data", inName);
    bwgCreate(sectionList, chromSizeHash, 256, 1024, TRUE, outName);
    lmCleanup(&lm);
}

struct hash *readCsizeHash(char *filename)
/* copied from lift.c */
{
    struct lineFile *lf = lineFileOpen(filename, TRUE);
    struct hash *cHash = hashNew(10);
    char *words[2];
    while (lineFileRowTab(lf, words))
	hashAddInt(cHash, words[0], sqlSigned(words[1]));
    lineFileClose(&lf);
    return cHash;
}

int main(int argc, char *argv[])
/* Process command line. */
{
    if (argc != 4)
	errAbort("bad running of bwmake");
    struct hash *sizeHash = readCsizeHash(argv[1]);
    writeBw(argv[2], argv[3], sizeHash);
    hashFree(&sizeHash);
    return 0;
}
