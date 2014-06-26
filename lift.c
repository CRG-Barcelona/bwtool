/* bwtool_lift - lift bigWig file */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <jkweb/common.h>
#include <jkweb/obscure.h>
#include <jkweb/linefile.h>
#include <jkweb/hash.h>
#include <jkweb/options.h>
#include <jkweb/sqlNum.h>
#include <jkweb/basicBed.h>
#include <jkweb/chain.h>
#include <jkweb/binRange.h>
#include <jkweb/bigWig.h>
#include <beato/bigs.h>
#include "bwtool.h"
#include "bwtool_shared.h"

#define NANUM sqrt(-1)

struct liftOverChromMap
/* Remapping information for one (old) chromosome */
{
    char *name;                 /* Chromosome name. */
    struct binKeeper *bk;       /* Keyed by old position, values are chains. */
};

void usage_lift()
/* Explain usage of the lift program */
{
errAbort(
  "bwtool lift - project data base-by-base into a new assembly using\n"
  "   a liftOver chain file from UCSC.\n"
  "usage:\n"
  "   bwtool lift old.bw[:chr:start-end] oldToNew.liftOver.chain.gz new.bw\n"
  "options:\n"
  "   -sizes=new.sizes       if set use the chrom.sizes file specified instead of\n"
  "                          gathering the size information from the chain.\n"
  "                          This is one way to restrict the chromosomes lifted to\n"
  "                          in the output.\n"
  "   -unlifted=file.bed     save all the regions from the input not lifted\n"
  );
}

struct hash *qSizeHash(char *chainfile)
/* read the chain file and figure out what the chromosome sizes are on the query end */
{
    struct lineFile *lf = lineFileOpen(chainfile, TRUE);
    struct chain *ch;
    struct hash *csizes = hashNew(10);
    while ((ch = chainRead(lf)) != NULL)
    {
	char *chrom = ch->qName;
	int size = ch->qSize;
	if (!hashLookup(csizes, chrom))
	    hashAddInt(csizes, chrom, size);
	chainFree(&ch);
    }
    lineFileClose(&lf);
    return csizes;
}

struct hash *readLiftOverMapChainHash(char *fileName)
/* taken from kent/src/hg/lib/liftOver.c */
/* Read map file into hashes. */
{
    struct hash *chainHash = hashNew(10);
    struct lineFile *lf = lineFileOpen(fileName, TRUE);
    struct chain *chain;
    struct liftOverChromMap *map;

    while ((chain = chainRead(lf)) != NULL)
    {
	if ((map = hashFindVal(chainHash, chain->tName)) == NULL)
	{
	    AllocVar(map);
	    map->bk = binKeeperNew(0, chain->tSize);
	    hashAddSaveName(chainHash, chain->tName, map, &map->name);
	}
	binKeeperAdd(map->bk, chain->tStart, chain->tEnd, chain);
    }
    lineFileClose(&lf);
    return chainHash;
}

void freeChainHashMap(void **pHashMapVal)
/* also frees the vals */
{
    struct liftOverChromMap **p = (struct liftOverChromMap **)pHashMapVal;
    struct liftOverChromMap *map = *p;
    binKeeperFree(&map->bk);
    freez(p);
}

struct hash *genomePbw(struct hash *qSizes)
/* make a parallel hash of pbws given the size hash, also keyed on chrom name */
{
    struct hash *pbwHash = newHash(10);
    struct hashEl *list = hashElListHash(qSizes);
    struct hashEl *el;
    const double na = NANUM;
    int i;
    for (el = list; el != NULL; el = el->next)
    {
	int size = ptToInt(el->val);
	struct perBaseWig *pbw = alloc_perBaseWig(el->name, 0, size);
	for (i = 0; i < pbw->len; i++)
	    pbw->data[i] = na;
	pbw->name = cloneString(el->name);
	pbw->strand[0] = '+';
	pbw->strand[1] = '\0';
	hashAdd(pbwHash, el->name, pbw);
    }
    hashElFreeList(&list);
    return pbwHash;
}

struct hash *readCsizeHash(char *filename)
/* read in a chrom sizes file */
{
    struct lineFile *lf = lineFileOpen(filename, TRUE);
    struct hash *cHash = hashNew(10);
    char *words[2];
    while (lineFileRowTab(lf, words))
	hashAddInt(cHash, words[0], sqlSigned(words[1]));
    lineFileClose(&lf);
    return cHash;
}

enum remapResult
{
    problem = 0,
    deleted = 1,
    duplicated = 2,
    lifted = 3
};

static int chainAliSize(struct chain *chain)
/* Return size of all blocks in chain. */
{
struct cBlock *b;
int total = 0;
for (b = chain->blockList; b != NULL; b = b->next)
    total += b->qEnd - b->qStart;
return total;
}

static boolean mapThroughChain(struct chain *chain, double minRatio,
	int *pStart, int *pEnd, struct chain **retSubChain,
	struct chain **retChainToFree)
/* Map interval from start to end from target to query side of chain.
 * Return FALSE if not possible, otherwise update *pStart, *pEnd. */
{
    struct chain *subChain = NULL;
    struct chain *freeChain = NULL;
    int s = *pStart, e = *pEnd;
    int oldSize = e - s;
    int newCover = 0;
    int ok = TRUE;

    chainSubsetOnT(chain, s, e, &subChain, &freeChain);
    if (subChain == NULL)
    {
	*retSubChain = NULL;
	*retChainToFree = NULL;
	return FALSE;
    }
    newCover = chainAliSize(subChain);
    if (newCover < oldSize * minRatio)
	ok = FALSE;
    else if (chain->qStrand == '+')
    {
	*pStart = subChain->qStart;
	*pEnd = subChain->qEnd;
    }
    else
    {
	*pStart = subChain->qSize - subChain->qEnd;
	*pEnd = subChain->qSize - subChain->qStart;
    }
    *retSubChain = subChain;
    *retChainToFree = freeChain;
    return ok;
}

enum remapResult remapBase(struct hash *chainHash, char *orig_chrom, int orig_base, char **dest_chrom, int *dest_base)
{
    struct liftOverChromMap *map = hashFindVal(chainHash, orig_chrom);
    struct binElement *list = NULL;
    struct chain *chainHit = NULL;
    struct chain *toFree;
    struct chain *subChain;
    int start = orig_base, end = start+1;
    if (map)
	list = binKeeperFind(map->bk, start, start+1);
    if (!list)
	return deleted;
    else if (list->next != NULL)
    {
	slFreeList(&list);
	return duplicated;
    }
    chainHit = list->val;
    if (!mapThroughChain(chainHit, 1, &start, &end, &subChain, &toFree))
    {
	slFreeList(&list);
	return problem;
    }
    chainFree(&toFree);
    *dest_base = start;
    *dest_chrom = chainHit->qName;
    slFreeList(&list);
    return lifted;
}

void do_pass1(struct metaBig *mb, struct hash *chainHash, struct hash *gpbw)
/* do the first pass in the destination pbws */
/* remap all the origen bases and increment from zero in */
/* the destination the number of times the base in the */
/* destination is mapped to. For the second pass we'll */
/* only use the ones that are 1 here.  Ones that are > 1 */
/* will be considered places where the destination is */
/* repeated */
{
    struct bed *section;
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbwList = perBaseWigLoadContinue(mb, section->chrom, section->chromStart, section->chromEnd);
	struct perBaseWig *pbw;
	for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
	{
	    int i;
	    for (i = 0; i < pbw->len; i++)
	    {
		char *dest_chrom = NULL;
		int dest_start = 0;
		enum remapResult rmr = remapBase(chainHash, pbw->chrom, pbw->chromStart + i, &dest_chrom, &dest_start);
		if (rmr == lifted)
		{
		    struct perBaseWig *dest_chrom_pbw = (struct perBaseWig *)hashFindVal(gpbw, dest_chrom);
		    if (dest_chrom_pbw)
		    {
			if (isnan(dest_chrom_pbw->data[dest_start]))
			    dest_chrom_pbw->data[dest_start] = 1.0;
			else
			    dest_chrom_pbw->data[dest_start] += 1.0;
		    }
		}
	    }
	}
	perBaseWigFreeList(&pbwList);
    }
}

void do_pass2(struct metaBig *mb, struct hash *chainHash, struct hash *gpbw)
/* so now that everything is either zero, 1.0, or more, make everything that isn't 1.0 an NA */
{
    const double na = NANUM;
    struct bed *section;
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbwList = perBaseWigLoadContinue(mb, section->chrom, section->chromStart, section->chromEnd);
	struct perBaseWig *pbw;
	for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
	{
	    int i;
	    for (i = 0; i < pbw->len; i++)
	    {
		char *dest_chrom = NULL;
		int dest_start = 0;
		enum remapResult rmr = remapBase(chainHash, pbw->chrom, pbw->chromStart + i, &dest_chrom, &dest_start);
		if (rmr == lifted)
		{
		    struct perBaseWig *dest_chrom_pbw = (struct perBaseWig *)hashFindVal(gpbw, dest_chrom);
		    if (dest_chrom_pbw && ((int)dest_chrom_pbw->data[dest_start] != 1))
			dest_chrom_pbw->data[dest_start] = na;
		}
	    }
	}
	perBaseWigFreeList(&pbwList);
    }
}

void do_final_pass(struct metaBig *mb, struct hash *chainHash, struct hash *gpbw, char *bad_file)
/* now everything is 1.0 or NA.  copy data into destination */
{
    struct bed *section;
    FILE *bad = NULL;
    if (bad_file)
	bad = mustOpen(bad_file, "w");
    for (section = mb->sections; section != NULL; section = section->next)
    {
	struct perBaseWig *pbwList = perBaseWigLoadContinue(mb, section->chrom, section->chromStart, section->chromEnd);
	struct perBaseWig *pbw;
	for (pbw = pbwList; pbw != NULL; pbw = pbw->next)
	{
	    int i;
	    for (i = 0; i < pbw->len; i++)
	    {
		char *dest_chrom = NULL;
		int dest_start = 0;
		enum remapResult rmr = remapBase(chainHash, pbw->chrom, pbw->chromStart + i, &dest_chrom, &dest_start);
		if (rmr == lifted)
		{
		    struct perBaseWig *dest_chrom_pbw = (struct perBaseWig *)hashFindVal(gpbw, dest_chrom);
		    if (dest_chrom_pbw && (!isnan(dest_chrom_pbw->data[dest_start])))
			dest_chrom_pbw->data[dest_start] = pbw->data[i];
		    else if (bad && isnan(dest_chrom_pbw->data[dest_start]))
			fprintf(bad, "%s\t%d\tmulti_mapped_%s_%d\n", pbw->chrom, pbw->chromStart+i, dest_chrom, dest_start);
		}
		else if (bad)
		{
		    if (rmr == duplicated)
			fprintf(bad, "%s\t%d\tduplicated_in_destination\n", pbw->chrom, pbw->chromStart+i);
		    else if (rmr == deleted)
			fprintf(bad, "%s\t%d\tdeleted_in_destination\n", pbw->chrom, pbw->chromStart+i);
		    else
			fprintf(bad, "%s\t%d\tproblem_lifting\n", pbw->chrom, pbw->chromStart+i);
		}
	    }
	}
	perBaseWigFreeList(&pbwList);
    }
    if (bad_file)
	carefulClose(&bad);
}

int pbwHashElCmp(const void *va, const void *vb)
{
    const struct perBaseWig *pbwa = (struct perBaseWig *)(*((struct hashEl **)va))->val;
    const struct perBaseWig *pbwb = (struct perBaseWig *)(*((struct hashEl **)vb))->val;
    return pbwb->len - pbwa->len;
}

void bwtool_lift(struct hash *options, char *favorites, char *regions, unsigned decimals,
		 enum wigOutType wot, char *bigfile, char *tmp_dir, char *chainfile, char *outputfile)
/* bwtool_lift - main for lifting program */
{
    struct hash *sizeHash = NULL;
    struct hash *chainHash = readLiftOverMapChainHash(chainfile);
    struct hash *gpbw = NULL;
    char *size_file = hashFindVal(options, "sizes");
    char *bad_file = hashFindVal(options, "unlifted");
    if (size_file)
	sizeHash = readCsizeHash(size_file);
    else
	sizeHash = qSizeHash(chainfile);
    gpbw = genomePbw(sizeHash);
    struct metaBig *mb = metaBigOpen_check(bigfile, tmp_dir, regions);
    char wigfile[512];
    safef(wigfile, sizeof(wigfile), "%s.tmp.wig", outputfile);
    FILE *out = mustOpen(wigfile, "w");
    struct hashEl *elList = hashElListHash(gpbw);
    struct hashEl *el;
    verbose(2,"starting first pass\n");
    do_pass1(mb, chainHash, gpbw);
    verbose(2, "starting second pass\n");
    do_pass2(mb, chainHash, gpbw);
    verbose(2,"starting final pass\n");
    do_final_pass(mb, chainHash, gpbw, bad_file);
    slSort(&elList, pbwHashElCmp);
    for (el = elList; el != NULL; el = el->next)
    {
	struct perBaseWig *pbw = (struct perBaseWig *)el->val;
	perBaseWigOutputNASkip(pbw, out, wot, decimals, NULL, FALSE, FALSE);
    }
    hashElFreeList(&elList);
    carefulClose(&out);
    hashFreeWithVals(&chainHash, freeChainHashMap);
    hashFreeWithVals(&gpbw, perBaseWigFree);
    writeBw(wigfile, outputfile, sizeHash);
    hashFree(&sizeHash);
    remove(wigfile);
    metaBigClose(&mb);
}
