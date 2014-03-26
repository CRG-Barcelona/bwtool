#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "common.h"
#include "obscure.h"
#include "hash.h"
#include "linefile.h"
#include "localmem.h"
#include "sqlNum.h"
#include "sig.h"
#include "basicBed.h"
#include "bigBed.h"
#include "bigWig.h"
#include "udc.h"
#include "rangeTree.h"
#include "metaBig.h"
#include "bigs.h"

#ifdef USE_BAM
#include "metaBigBam.h"
#endif

/********** functions that may need going into basicBed.ch some day ***************/

struct bed6 *readBed6(char *file)
/* read from a file */
{
    char *words[6];
    if (!fileExists(file))
	errAbort("Can't find file: %s", file);
    struct lineFile *lf = lineFileOpen(file, TRUE);
    struct bed6 *list = NULL;
    int num_words = 0;
    while ((num_words = lineFileChopTab(lf, words)) > 0)
    {
	if (num_words < 6)
	    errAbort("Expecting BED-6 formatted file (%s) but there are only %d columns on line %d", file, num_words, lf->lineIx);
	struct bed6 *newb;
	AllocVar(newb);
	newb->chrom = cloneString(words[0]);
	newb->chromStart = sqlSigned(words[1]);
	newb->chromEnd = sqlSigned(words[2]);
	newb->name = cloneString(words[3]);
	newb->score = sqlSigned(words[4]);
	newb->strand[0] = words[5][0];
	newb->strand[1] = '\0';
	slAddHead(&list, newb);
    }
    slReverse(&list);
    lineFileClose(&lf);
    return list;
}

struct bed6 *readBed6SoftAndSize(char *file, int *orig_size)
/* read from a file.  If it's missing fields, fill the bed */
{
    char *words[6];
    if (!fileExists(file))
	errAbort("Can't find file: %s", file);
    struct lineFile *lf = lineFileOpen(file, TRUE);
    struct bed6 *list = NULL;
    int num_words = 0;
    int size = 0;
    while ((num_words = lineFileChopTab(lf, words)) > 0)
    {
	if (num_words < 3)
	    errAbort("Expecting BED-3 formatted file (%s) but there are only %d columns on line %d", file, num_words, lf->lineIx);
	struct bed6 *newb;
	AllocVar(newb);
	if (size < num_words)
	    size = num_words;
	newb->chrom = cloneString(words[0]);
	newb->chromStart = sqlSigned(words[1]);
	newb->chromEnd = sqlSigned(words[2]);
	if (num_words < 4)
	    newb->name = cloneString(".");
	else
	    newb->name = cloneString(words[3]);
	if (num_words < 5)
	    newb->score = 0;
	else
	    newb->score = sqlSigned(words[4]);
	if (num_words < 6)
	    newb->strand[0] = '+';
	else
	    newb->strand[0] = words[5][0];
	newb->strand[1] = '\0';
	slAddHead(&list, newb);
    }
    slReverse(&list);
    lineFileClose(&lf);
    if (orig_size)
	*orig_size = size;
    return list;
}

struct bed6 *readBed6Soft(char *file)
/* read from a file.  If it's missing fields, fill the bed */
{
    return readBed6SoftAndSize(file, NULL);
}

static struct bed *readAtLeastBed3(char *file)
/* read from a file */
{
    char *words[3];
    if (!fileExists(file))
	errAbort("Can't find file: %s", file);
    struct lineFile *lf = lineFileOpen(file, TRUE);
    struct bed *list = NULL;
    int num_words = 0;
    while ((num_words = lineFileChopTab(lf, words)) > 0)
    {
	if (num_words < 3)
	    errAbort("Expecting BED-3 formatted file (%s) but there are only %d columns on line %d", file, num_words, lf->lineIx);
	struct bed *newb;
	AllocVar(newb);
	newb->chrom = cloneString(words[0]);
	newb->chromStart = sqlSigned(words[1]);
	newb->chromEnd = sqlSigned(words[2]);
	slAddHead(&list, newb);
    }
    slReverse(&list);
    lineFileClose(&lf);
    return list;
}

void bed6Free(struct bed6 **pEl)
/* Free a single dynamically allocated bed such as created
 * with bedLoad(). */
{
struct bed6 *el;
if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freez(pEl);
}

void bed6FreeList(struct bed6 **pList)
/* Free a list of dynamically allocated bed's */
{
struct bed6 *el, *next;
for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    bed6Free(&el);
    }
*pList = NULL;
}

void pairbedFree(struct pairbed **pEl)
/* free up a single pairbed */
{
    struct pairbed *el;
    if ((el = *pEl) == NULL) return;
    freeMem(el->chrom);
    freeMem(el->mChrom);
    freeMem(el->name);
    freez(pEl);
}

void pairbedFreeList(struct pairbed **pList)
/* free up a list of them */
{
    struct pairbed *el, *next;
    for (el = *pList; el != NULL; el = next)
    {
	next = el->next;
	pairbedFree(&el);
    }
    *pList = NULL;    
}

enum metaBigFileType isBigWigOrBed(char *filename)
/* Peak at a file to see if it's bigWig */
{
    enum metaBigFileType ret = isNotBig;
    bits32 magic;
    struct udcFile *udc;
    udc = udcFileOpen(filename, udcDefaultDir());
    /* Read magic number at head of file and use it to see if we are proper file type, and
     * see if we are byte-swapped. */
    udcMustRead(udc, &magic, sizeof(magic));
    if (magic == bigWigSig)
	ret = isaBigWig;
    else if (magic == bigBedSig)
	ret = isaBigBed;
    if (ret != isNotBig)
    {
	udcFileClose(&udc);
	return ret;
    }
    magic = byteSwap32(magic);
    if (magic == bigWigSig)
	ret = isaBigWig;
    else if (magic == bigBedSig)
	ret = isaBigBed;
    udcFileClose(&udc);
    return ret;
}

enum metaBigFileType sniffBigFile(char *filename)
/* try to figure out what type of big file it is. */
{
    enum metaBigFileType ft = isBigWigOrBed(filename);
    if ((ft == isaBigBed) || (ft == isaBigWig))
	return ft;
#ifdef USE_BAM
    if (isBamWithIndex(filename))
	return isaBam;
#endif
    return isNotBig;
}

/**** chrom-sizes hash code *****/ 

static struct hash *bbiChromSizes(struct bbiFile *bbi)
/* return the hash of chrom sizes from the bigBed/bigWig */
{
    struct bbiChromInfo *cList = bbiChromList(bbi);
    struct bbiChromInfo *c;
    struct hash *cHash = newHash(10);
    for (c = cList; c != NULL; c = c->next)
	hashAddInt(cHash, c->name, (int)c->size);
    bbiChromInfoFreeList(&cList);
    return cHash;
}


static struct bed *sectionsFromChromSizes(struct hash *chromSizes)
/* return a list of bed3s used as sections using the whole chromosome sizes */
{
    struct bed *bedList = NULL;
    struct hashEl *list = hashElListHash(chromSizes);
    struct hashEl *el;
    for (el = list; el != NULL; el = el->next)
    {
	struct bed *bed;
	AllocVar(bed);
	bed->chrom = cloneString(el->name);
	bed->chromStart = 0;
	bed->chromEnd = (unsigned)ptToInt(el->val);
	slAddHead(&bedList, bed);
    }
    slReverse(&bedList);
    hashElFreeList(&list);
    return bedList;
}

/****** functions to deal with the filename *****************/

static boolean parseMetaBigFileName(char *original, char **pRemotePart, char **pFullFilePart, char **pJustFilePart, char **pSectionPart)
/* split the filename string into up to three pieces: remote/file/sections parts */
{
    boolean remote = FALSE;
    char *baseFileName = NULL;
    char *colonInFileName = NULL;
    if (original == NULL)
	errAbort("NULL filename error, something went wrong");
    if (strstr(original, "tp://") || strstr(original, "https://"))
	remote = TRUE;
    /* find right-most '/' */
    baseFileName = strrchr(original, '/');
    if (baseFileName == NULL)
    {
	baseFileName = original;
	*pRemotePart = NULL;
    }
    else if (remote && ((baseFileName - original == 7) || (baseFileName - original == 6)))
	errAbort("malformed URL");
    else if (remote)
    {
	*pRemotePart = cloneStringZ(original, baseFileName+1 - original);
	baseFileName += 1;
    }
    else
	baseFileName += 1;
    colonInFileName = strchr(baseFileName, ':');
    if (colonInFileName != NULL)    
    /* this means there is a section portion of the filename */
    {
	*pFullFilePart = cloneStringZ(original, colonInFileName - original);
	*pJustFilePart = cloneStringZ(baseFileName, colonInFileName - baseFileName);
	if (colonInFileName+1 == NULL)
	    errAbort("Sectioning operator \':\' used after filename but left empty"); 
	*pSectionPart = cloneString(colonInFileName+1);
    }
    else
    {
	*pFullFilePart = cloneString(original);
	*pJustFilePart = cloneString(baseFileName);
	*pSectionPart = NULL;
    }
    return remote;
}

static struct bed *sectionToBed(char *sectionString, struct hash *chromSizes)
/* deconstruct the "chr:start-end" into what's needed for the */
/* bed.  Also allow just "chr" or "chr:start+size" */
{
    struct bed *sec = NULL;
    char *colon = NULL;
    AllocVar(sec);
    colon = strchr(sectionString, ':');
    if (colon == NULL)
    /* only a chrom */
    {
	sec->chrom = cloneString(sectionString);
	sec->chromStart = 0;
	sec->chromEnd = (unsigned)hashIntVal(chromSizes, sec->chrom);
    }
    else
    {
	unsigned chromSize = 0;
	char *minusOrPlus = NULL;
	char *start = colon+1;
	char *end = NULL;
	if (start == NULL)
	    errAbort("malformed range (nothing after ':')");
	sec->chrom = cloneStringZ(sectionString, colon - sectionString);
	chromSize = (unsigned)hashIntVal(chromSizes, sec->chrom);
	minusOrPlus = strchr(start, '-');
	if (minusOrPlus != NULL)
	{
	    if (minusOrPlus == start)
		sec->chromStart = 0;
	    else
	    {
		*minusOrPlus = '\0';
		sec->chromStart = sqlUnsigned(start);
		*minusOrPlus = '-';
	    }
	    end = minusOrPlus + 1;
	    if ((end != NULL) && (strlen(end) > 0))
		sec->chromEnd = sqlUnsigned(end);
	    else
		sec->chromEnd = chromSize;
	}
	else
	{
	    /* check for '+' */
	    minusOrPlus = strchr(start, '+');
	    if (minusOrPlus != NULL)
	    {
		if (minusOrPlus == start)
		    sec->chromStart = 0;
		else
		{
		    *minusOrPlus = '\0';
		    sec->chromStart = sqlUnsigned(start);
		    *minusOrPlus = '+';
		}
		end = minusOrPlus + 1;
		if (end != NULL)
		    sec->chromEnd = sqlUnsigned(end) + sec->chromStart;
	    }
	    else 
		errAbort("malformed range %s", sectionString);
	}
	if (sec->chromStart >= sec->chromEnd)
	    errAbort("bad range specified: start >= end with %s", sectionString);
	if (sec->chromEnd > chromSize)
	    errAbort("bad range specified: end > %s size (%u) for %s", sec->chrom, chromSize, sectionString);
    }
    return sec;
}

static struct bed *parseSectionString(char *sectionString, struct hash *chromSizes)
/* returns a list of sections represented as beds as anything */
/* past a ':' given in the file name.  given string should only be that portion.*/
/* passing a NULL sectionString produces the entire chromosomes from the */
/* chromSizes hash table. */
{
    char *sections[1024];
    struct bed *bedList = NULL;
    if (sectionString)
    {
	int numSections = 0;
	int i;
	char *s = cloneString(sectionString);
	stripChar(s, ',');
	numSections = chopString(s, ";", sections, ArraySize(sections));
	for (i = 0; i < numSections; i++)
	{
	    struct bed *sec = sectionToBed(sections[i], chromSizes);
	    slAddHead(&bedList, sec);
	}
	slReverse(&bedList);
    }
    else
	bedList = sectionsFromChromSizes(chromSizes);
    return bedList;
}

/** General metaBig functions **/


/************** BAM-specific retrival code (could maybe go into bamFile.ch some day) *************************/


/************** bigBed fetching code ******************************/

static struct bed6 *bed6FromBigBedInterval(struct bigBedInterval *bbi, struct metaBigBed6Helper helper) 
/* convert a bigBedInterval into something a bit more userful */
{
    struct bed6 *bed;
    struct lm *lm = helper.lm;
    struct metaBig *mb = helper.mb;
    int start = bbi->start;
    int end = bbi->end;
    char *words[32];
    char sBuf[16], eBuf[16];
    int numFields = bigBedIntervalToRow(bbi, helper.chrom, sBuf, eBuf, words, 32);
    char strand = 0;
    if (numFields < 6)
	errAbort("Expecting six fields in the bigBed and there are %d", numFields);
    /* decide whether to add this or not now, because if so, more lm memory is to be alloc'd */
    strand = words[5][0];
    if ((mb->strand != 0) && (mb->strand != strand))
    	return NULL;
    if (strand == '+')
    {
	start += mb->shift;
	if (mb->length > 0)
	    end = start + mb->length;
	start += mb->shift;
	end += mb->shift;
    }
    else
    {
	if (mb->length > 0)
	    start = end - mb->length;;
	start -= mb->shift;
	end -= mb->shift;
    }
    /* filter out ones that go out-of-bounds after shifting/extending */
    if ((start >= end) || (start < 0) || (end > helper.chromSize))
    	return NULL;
    lmAllocVar(lm, bed);
    bed->chrom = helper.chrom;
    bed->chromStart = start;
    bed->chromEnd = end;
    if (mb->nameType == basicName)
	bed->name = lmCloneString(lm, words[3]);
    else 
	bed->name = helper.dot;
    bed->score = sqlUnsigned(words[4]);
    bed->strand[0] = strand;
    return bed;
}

static struct bed6 *bigBedBed6Fetch(struct metaBig *mb, char *chrom, unsigned start, unsigned end, struct lm *lm)
/* the main fetcher */
{
    struct metaBigBed6Helper helper;
    struct bigBedInterval *ints = bigBedIntervalQuery(mb->big.bbi, chrom, (bits32)start, (bits32)end, 0, lm);
    struct bigBedInterval *ival;
    helper.lm = lm;
    helper.chrom = chrom;
    helper.chromSize = hashIntVal(mb->chromSizeHash, chrom);
    helper.dot = ".";
    helper.bedList = NULL;
    helper.mb = mb;
    for (ival = ints; ival != NULL; ival = ival->next)
    {
	struct bed6 *newbed = bed6FromBigBedInterval(ival, helper);
	if (newbed)
	    slAddHead(&helper.bedList, newbed);
    }
    slReverse(&helper.bedList);
    return helper.bedList;
}

static long bigBedCount(struct metaBig *mb, char *chrom, unsigned start, unsigned end, boolean fifty)
/* main counter for bigBed.  I wonder if this can be done without fetching items */
{
    struct lm *lm = lmInit(0);
    struct bigBedInterval *ints = bigBedIntervalQuery(mb->big.bbi, chrom, (bits32)start, (bits32)end, 0, lm);
    if (fifty)
    {
	struct bigBedInterval *new_list = NULL;
	struct bigBedInterval *one = NULL;
	while ((one = slPopHead(&ints)) != NULL)
	{
	    unsigned middle = (unsigned)(((double)one->end-one->start)/2)+one->start;
	    if ((middle >= start) && (middle < end))
		slAddHead(&new_list, one);
	}
	ints = new_list;
    }
    slReverse(&ints);
    int count = slCount(ints);
    lmCleanup(&lm);
    return (long)count;
}

static struct bed *subset_beds(char *sectionString, struct bed **pRegions, struct hash *chromHash)
/* in the situation where both a regions bed file is given AND the filename specifies subsections, */
/* intersect the two.  For simplictity sake,  */
{
    struct bed *fname_ranges = parseSectionString(sectionString, chromHash);
    struct bed *bed;
    struct bed *subset = NULL;
    struct bed *regions = *pRegions;
    slSort(&fname_ranges, bedCmp);
    bed = fname_ranges;
    while (bed != NULL)
    {
	/* each iteration of the loop should be a separate chrom */
	struct bed *region;
	struct rbTree *tree = rangeTreeNew();
	while ((bed != NULL) && (bed->next != NULL) && (sameString(bed->chrom, bed->next->chrom)))
	{
	    rangeTreeAdd(tree, bed->chromStart, bed->chromEnd);
	    bed = bed->next;
	}
	rangeTreeAdd(tree, bed->chromStart, bed->chromEnd);
	/* now we're at a point that we're dealing only with one chromosome. */
	for (region = regions; region != NULL; region = region->next)
	{
	    if (sameString(region->chrom, bed->chrom) && rangeTreeOverlaps(tree, region->chromStart, region->chromEnd)
		&& rangeTreeFindEnclosing(tree, region->chromStart, region->chromEnd))
	    {
		struct bed *clone = cloneBed(region);
		slAddHead(&subset, clone);
	    }
	    else if (sameString(region->chrom, bed->chrom) && rangeTreeOverlaps(tree, region->chromStart, region->chromEnd))
		    errAbort("range specified in file overlaps but is not contained by range specified on command-line");
	}
	rangeTreeFree(&tree);
	bed = bed->next;
    }
    if (subset == NULL)
    {
	errAbort("no ranges specified in file were contained in ranges specified on command-line");
    }
    slReverse(&subset);
    bedFreeList(&fname_ranges);
    bedFreeList(pRegions);
    return subset;
}

static struct bed *regionsLoad(char *sectionsBed)
/* return a bed3 list of regions for times when -regions is used. */
/* If the filename has a comma then a number, then take just that line */
{
    struct bed *list = NULL;
    unsigned ix = 0;
    char *comma = NULL;
    if (strchr(sectionsBed, ','))
    {
	char *number_part = chopPrefixAt(sectionsBed, ',');
	if (number_part)
	    ix = sqlUnsigned(number_part);
    }
    list = readAtLeastBed3(sectionsBed);
    if (list && (ix > 0))
    {
	struct bed *single = slElementFromIx(list, ix-1);
	if (single)
	{
	    struct bed *rem;
	    while ((rem = slPopHead(&list)) != single)
		bedFree(&rem);
	    rem = single->next;
	    bedFreeList(&rem);
	    single->next = NULL;
	    list = single;
	}
    }
    return list;
}

/************** the "public" functions from the .h  ***************/

struct metaBig *metaBigOpen(char *fileOrUrlwSections, char *sectionsBed)
/* load a file or URL with or without sectioning */
/* if it's a bam, load the index. */
{
    struct metaBig *mb;
    char *fullFileName = NULL;
    char *remoteDir = NULL;
    char *baseFileName = NULL;
    char *sections = NULL;
    AllocVar(mb);
    mb->originalFileName = cloneString(fileOrUrlwSections);
    /* first deal with filename and separate URL/file/sections */
    mb->isRemote = parseMetaBigFileName(fileOrUrlwSections, &remoteDir, &fullFileName, &baseFileName, &sections);
    mb->fileName = fullFileName;
    mb->baseFileName = baseFileName;
    mb->remoteSiteAndDir = remoteDir;
    /* sniff the file */
    mb->type = sniffBigFile(mb->fileName);
    /* depending on the type, open the files and get the chrom-size hash different ways */
    if (mb->type == isaBigBed)
    {
	mb->big.bbi = bigBedFileOpen(mb->fileName);
	mb->chromSizeHash = bbiChromSizes(mb->big.bbi);
	mb->numReads = bigBedItemCount(mb->big.bbi);
    }
#ifdef USE_BAM
    else if (mb->type == isaBam)
    {
	char *anotherFileName = NULL;
	mb->chromSizeHash = bamChromSizes(mb->fileName);
	mb->header = bamGetHeaderOnly(mb->fileName);
	mb->big.bam = bamOpen(mb->fileName, &anotherFileName);
	if (!sameString(mb->fileName, anotherFileName))
	    mb->fileName = anotherFileName;
	/* Also need to load the index since it's a bam */
	mb->idx = bam_index_load(mb->fileName);
	metaBigBamFlagCountsInit(mb);
    }
#endif
    else if (mb->type == isaBigWig)
    {
	mb->big.bbi = bigWigFileOpen(mb->fileName);
	mb->chromSizeHash = bbiChromSizes(mb->big.bbi);	
    }
    else 
    {
	/* maybe I should free some stuff up here */
	if (fullFileName)
	    freeMem(fullFileName);
	if (remoteDir)
	    freeMem(remoteDir);
	if (baseFileName)
	    freeMem(baseFileName);
	if (sections)
	    freeMem(sections);
	freez(&mb);
	return NULL;
    }
    if (sectionsBed && sections)
    {
	struct bed *regions = (fileExists(sectionsBed)) ? regionsLoad(sectionsBed) : parseSectionString(sectionsBed, mb->chromSizeHash);
	struct bed *subsets = subset_beds(sections, &regions, mb->chromSizeHash);
	mb->sections = subsets;
    }
    else if (sectionsBed)
    {
	mb->sections = (fileExists(sectionsBed)) ? regionsLoad(sectionsBed) : parseSectionString(sectionsBed, mb->chromSizeHash);
    }
    else
	mb->sections = parseSectionString(sections, mb->chromSizeHash);
    return mb;
}

void metaBigClose(struct metaBig **pMb)
/* close the file and free up everything. */
{
    struct metaBig *mb = *pMb;
    hashFree(&mb->chromSizeHash);
    if (mb->rgList)
	hashFree(&mb->rgList);
    if (mb->sections)
	bedFreeList(&mb->sections);
    if (mb->originalFileName)
	freeMem(mb->originalFileName);
    if (mb->fileName)
	freeMem(mb->fileName);
    if (mb->baseFileName)
	freeMem(mb->baseFileName);
    if (mb->remoteSiteAndDir)
	freeMem(mb->remoteSiteAndDir);
#ifdef USE_BAM
    if (mb->idx)
	bam_index_destroy(mb->idx);
#endif
    if (mb->type == isaBigBed)
	bigBedFileClose(&mb->big.bbi);
#ifdef USE_BAM
    else if (mb->type == isaBam)
	samclose(mb->big.bam);
#endif
    else
	bigWigFileClose(&mb->big.bbi);
#ifdef USE_BAM
    if (mb->header)
	bam_header_destroy(mb->header);
#endif
    freez(pMb);
}

void metaBigSetPositionalOptions(struct metaBig *mb, int length, int shift, char strand)
/* set retrieval options for the metaBig. */
{
    mb->length = length;
    mb->shift = shift;
    mb->strand = strand;
}

void metaBigSetNameOption(struct metaBig *mb, enum metaBigNameType nameType)
/* set the naming option for retrieval */
{
    mb->nameType = nameType;
}

void metaBigToggleBadRegions(struct metaBig *mb, boolean useBad)
/* setting this to TRUE permits the retrieval of reads marked as from "bad regions" */
/* from the metaBig (only relevant with BAM files using the ZR tag) */
{
    mb->includeBadRegions = TRUE;
}

void metaBigToggleBRead(struct metaBig *mb, boolean useB)
/* setting this to TRUE permits the retrieval of reads marked as B-reads from the metaBig */
/* (only relevant with BAM files using the ZL tag) */
{
    mb->includeB = TRUE;
}

void metaBigSetRgList(struct metaBig *mb, char *commaSep, boolean isBlackList)
/* sets the whitelist (or blacklist) */
{
    struct slName *commas = slNameListFromComma(commaSep);
    struct slName *rgName; 
    struct hash *rgHash = newHash(10);
    for (rgName = commas; rgName != NULL; rgName = rgName->next)
	hashAdd(rgHash, rgName->name, rgName->name);
    slNameFreeList(&commas);
    mb->rgList = rgHash;
    mb->rgListIsBlack = isBlackList;
}

struct bed6 *metaBigBed6Fetch(struct metaBig *mb, char *chrom, unsigned start, unsigned end, struct lm *lm)
/* the main fetcher */
{
    if (mb->type == isaBigBed)
	return bigBedBed6Fetch(mb, chrom, start, end, lm);
#ifdef USE_BAM
    if (mb->type == isaBam)
	return bamBed6Fetch(mb, chrom, start, end, lm);
#endif
    return NULL;
}

struct pairbed *metaBigPairbedFetch(struct metaBig *mb, char *chrom, unsigned start, unsigned end, struct lm *lm)
/* fetcher for pairbed for Hi-C */
{
    if (mb->type == isaBigBed)
	errAbort("bigBed not supported yet");
#ifdef USE_BAM
    if (mb->type == isaBam)
	return bamPairbedFetch(mb, chrom, start, end, lm);
#endif
    return NULL;
}

long metaBigCount(struct metaBig *mb, char *chrom, unsigned start, unsigned end)
/* the main counter */
{
    if (mb->type == isaBigBed)
	return bigBedCount(mb, chrom, start, end, FALSE);
#ifdef USE_BAM
    if (mb->type == isaBam)
	return bamCount(mb, chrom, start, end);
#endif
    return -1;
}

long metaBigFiftyCount(struct metaBig *mb, char *chrom, unsigned start, unsigned end)
/* the main counter */
{
    if (mb->type == isaBigBed)
	return bigBedCount(mb, chrom, start, end, TRUE);
#ifdef USE_BAM
    if (mb->type == isaBam)
	errAbort("-fifty option not possible with bams right now");
#endif
    return -1;
}

long metaBigNumItems(struct metaBig *mb, boolean verbose)
/* return the total number of items in a bigBed or BAM */
/* used on a bigWig will return 0 */
/* unfortunately this is a loop through the entire file basically. */
/* nicer would be something that just glances at the index, but doing that */
/* might count items that would be filtered out upon fetching. */
{
    long sum = 0;
    struct bed *section;
    struct bed *chroms = NULL;
    if (mb->type == isaBigWig)
	return 0;
    else if (mb->type == isaBigBed)
	return (long)bigBedItemCount(mb->big.bbi);
    else
	chroms = sectionsFromChromSizes(mb->chromSizeHash);
    for (section = chroms; section != NULL; section = section->next)
    {
	struct lm *lm = lmInit(0);
	struct bed6 *list = metaBigBed6Fetch(mb, section->chrom, section->chromStart, section->chromEnd, lm);
	int num = slCount(list);
	if (verbose)
	    printf("Number of items in %s of %s: %d\n", section->chrom, mb->fileName, num);
	sum += num;
	lmCleanup(&lm);
    }
    bedFreeList(&chroms);
    return sum;
}
