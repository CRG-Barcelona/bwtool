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
#include "bamFile.h"
#include "rangeTree.h"
#include "metaBig.h"
#include "bigs.h"

#include "sam.h"
#include "bam.h"

/********** functions that may need going into basicBed.ch some day ***************/

struct bed6 *readBed6(char *file)
/* read from a file */
{
    char *words[6];
    struct lineFile *lf = lineFileOpen(file, TRUE);
    struct bed6 *list = NULL;
    while (lineFileRowTab(lf, words))
    {
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

/**** file-sniffing code ****/

typedef struct {  
    int beg, end;  
    samfile_t *in;  
} tmpstruct_t;

static boolean isBamWithIndex(char *file)
/* attempt to open a BAM file and its index. if everything is ok, */
/* return TRUE before closing. */
{
    tmpstruct_t tmp;
    bam_index_t *idx;
    tmp.beg = 0; 
    tmp.end = 0x7fffffff;  
    tmp.in = samopen(file, "rb", 0);  
    if (!tmp.in)
	return FALSE;
    idx = bam_index_load(file);
    if (!idx)
	return FALSE;
    bam_index_destroy(idx);
    samclose(tmp.in);
    return TRUE;
}

static enum metaBigFileType isBigWigOrBed(char *filename)
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
    magic = byteSwap32(magic);
    if (magic == bigWigSig)
	ret = isaBigWig;
    else if (magic == bigBedSig)
	ret = isaBigBed;
    udcFileClose(&udc);
    return ret;
}

static enum metaBigFileType sniffBigFile(char *filename)
/* try to figure out what type of big file it is. */
{
    enum metaBigFileType ft = isBigWigOrBed(filename);
    if ((ft == isaBigBed) || (ft == isaBigWig))
	return ft;
    if (isBamWithIndex(filename))
	return isaBam;
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

static bam_header_t *bamGetHeaderOnly(char *file)
/* just get the header, close the BAM */
{
    bamFile fp;
    bam_header_t *header;
    fp = bam_open(file, "r");
    if (!fp)
	errAbort("failed to open BAM file %s", file);
    header = bam_header_read(fp);
    bam_close(fp);
    return header;
}

static struct hash *bamChromSizes(char *bamfile)
/* return a list of chromosome sizes from a bam file */
{
    int i;
    struct hash *cHash = newHash(10);
    bam_header_t *header = bamGetHeaderOnly(bamfile);
    for (i = 0; i < header->n_targets; i++)
	hashAddInt(cHash, header->target_name[i], (int)header->target_len[i]);
    bam_header_destroy(header);
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
	numSections = chopString(s, ",", sections, ArraySize(sections));
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

struct metaBigCountHelper
/* Helper structure for count-only retrival */
{
    long count;
    int chromSize;
    struct metaBig *mb;
};

static int bamAddCount(const bam1_t *bam, void *data)
/* bam_fetch() calls this on each bam alignment retrieved.  Translate each bam
 * into a bed. */
{
    struct metaBigCountHelper *helper = (struct metaBigCountHelper *)data;
    struct metaBig *mb = helper->mb;
    const bam1_core_t *core = &bam->core;
    uint8_t *zr_tag = bam_aux_get(bam, "ZR");
    uint8_t *zl_tag = bam_aux_get(bam, "ZL");
    uint8_t *zd_tag = bam_aux_get(bam, "ZD");
    uint8_t *rg_tag = bam_aux_get(bam, "RG");
    char *rg = (rg_tag != NULL) ? bam_aux2Z(rg_tag) : NULL;
    /* first find out if we're to skip this read because there's a strand choice */
    /* in the options */
    char strand = '+';
    int start, end;
    /* don't consider if filtered */
    if ((core->flag & BAM_FQCFAIL) || (core->flag & BAM_FDUP))
	return 0;
    if (!mb->includeB && zl_tag)
	return 0;
    if (!mb->includeBadRegions && zr_tag)
	return 0;
    if ((mb->mapQ > 0) && ((int)core->qual < mb->mapQ))
	return 0;
    /* check whitelist/blacklist */
    if ((rg) && (mb->rgList) && ((!mb->rgListIsBlack && !hashLookup(mb->rgList, rg)) || (mb->rgListIsBlack && hashLookup(mb->rgList, rg))))
	return 0;
    if (bamIsRc(bam))
	strand = '-';
    if ((mb->useDupes > 1) && (zd_tag))
    {
	char *s = bam_aux2Z(zd_tag);
	char *num = chopPrefix(s);
	unsigned n = sqlUnsigned(num);
	freeMem(s);
	if (n > mb->useDupes)
	    return 0;
    }
    /* check single or paired-end */
    if ((!(core->flag&BAM_FPAIRED)) || (core->flag & BAM_FMUNMAP))
    /* the read isn't paired or its mate is unmapped */ 
    {
	if ((mb->strand != 0) && (mb->strand != strand))
	    return 0;
	start = core->pos;
	/* do shifts, extensions */
	if (strand == '+')
	{
	    if (mb->length > 0)
		end = start + mb->length;
	    else if (core->isize > 0)
		end = start + core->isize;
	    else 
		end = bam_calend(core, bam1_cigar(bam));
	}
	else
	{
	    end = bam_calend(core, bam1_cigar(bam));
	    if (mb->length > 0)
		start = end - mb->length;
	    else if (core->isize > 0)
		start = end - core->isize;
	}
	/* filter out ones that go out-of-bounds after shifting/extending */
	if ((start >= end) || (start < 0) || (end > helper->chromSize))
	    return 0;
	/* ok fine, we're adding it to the bed. */
	helper->count++;
    }
    else
    /* paired end */
    {
	if (core->flag & BAM_FPROPER_PAIR)
	{
	    /* only take one of the paired reads.  the other will be skipped */
	    if ((int)core->isize > 0)
	    {
		int mid;
		start = core->pos;
		end = start + core->isize;
		mid = core->pos + (core->isize / 2);
		if (mb->length > 0)
		{
		    start = mid - mb->length/2;
		    end = start + mb->length;
		}
		if ((start < 0) || (end > helper->chromSize))
		    return 0;
		helper->count++;
	    }
	}
    }
    return 0;
}

static long bamCount(struct metaBig *mb, char *chrom, unsigned start, unsigned end)
/* the main fetcher */
{
    struct metaBigCountHelper helper;
    char posForBam[256];
    helper.mb = mb;
    helper.count = 0;
    helper.chromSize = hashIntVal(mb->chromSizeHash, chrom);
    safef(posForBam, sizeof(posForBam), "%s:%d-%d", chrom, start, end);
    bamFetchAlreadyOpen(mb->big.bam, mb->idx, mb->fileName, posForBam, bamAddCount, &helper);
    return helper.count;
}

/************** BAM-specific retrival code (could maybe go into bamFile.ch some day) *************************/

struct metaBigBed6Helper
/* Helper structure to help get sam alignments out of a BAM file */
    {
    struct lm *lm;	/* Local memory pool to allocate into */
    char *chrom;	/* Chromosome name */
    char *dot;          /* When name='.', this remains dot. */
    int chromSize;      /* size of chromosome */
    struct metaBig *mb; /* the original metaBig with all the options and everything */
    struct bed6 *bedList;	/* List of alignments. */
    };

struct metaBigPairbedHelper
{
    struct lm *lm;
    char *chrom;
    char *dot;
    int chromSize;
    bam_header_t *header;
    struct metaBig *mb;
    struct pairbed *pbList;
};


static char *lmBamGetQuerySequence(const bam1_t *bam, boolean useStrand, struct lm *lm)
/* Allocate and return the nucleotide sequence encoded in bam.  The BAM format 
 * reverse-complements query sequence when the alignment is on the - strand,
 * so if useStrand is given we rev-comp it back to restore the original query 
 * sequence. */
{
const bam1_core_t *core = &bam->core;
int qLen = core->l_qseq;
char *qSeq;
lmAllocArray(lm, qSeq, qLen+1);
bamUnpackQuerySequence(bam, useStrand, qSeq);
return qSeq;
}

static char *lmBamGetQuality(const bam1_t *bam, boolean useStrand, struct lm *lm)
/* Get quality */
{
const bam1_core_t *core = &bam->core;
int qLen = core->l_qseq;
char *qSeq;
int i;
UBYTE *arr = bamGetQueryQuals(bam, useStrand);
lmAllocArray(lm, qSeq, qLen+1);
for (i = 0; i < qLen; i++)
    qSeq[i] = (char)(arr[i]+33);
qSeq[qLen] = '\0';
freeMem(arr);
return qSeq;
}

static char *lmBamGetOriginalName(const bam1_t *bam, struct lm *lm)
/* Allocate and return the nucleotide sequence encoded in bam.  The BAM format 
 * reverse-complements query sequence when the alignment is on the - strand,
 * so if useStrand is given we rev-comp it back to restore the original query 
 * sequence. */
{
    return lmCloneString(lm, bam1_qname(bam));
}

static boolean filterBam(struct metaBig *mb, const bam1_t *bam, const bam1_core_t *core)
/* most of the standard filtering stuff for each bam */
{
    uint8_t *rg_tag = bam_aux_get(bam, "RG");
    char *rg = (rg_tag != NULL) ? bam_aux2Z(rg_tag) : NULL;
    if ((mb->pe) && !(core->flag & BAM_FPROPER_PAIR))
	return TRUE;
    if ((mb->se) && (core->flag & BAM_FPROPER_PAIR))
	return TRUE;
    if ((core->flag & BAM_FQCFAIL) && (mb->includeB || mb->includeBadRegions || mb->includeMissingRestrSite || mb->useMapQ))
    {
	uint8_t *zr_tag = bam_aux_get(bam, "ZR");
	uint8_t *zl_tag = bam_aux_get(bam, "ZL");
	uint8_t *ze_tag = bam_aux_get(bam, "ZE");
	if (!mb->includeB && zl_tag)
	    return TRUE;
	if (!mb->includeBadRegions && zr_tag)
	    return TRUE;
	if (!mb->includeMissingRestrSite && ze_tag)
	    return TRUE;
	if ((mb->mapQ > 0) && ((int)core->qual < mb->mapQ))
	    return TRUE;
    }
    else if (core->flag & BAM_FQCFAIL)
	return TRUE;
    if ((core->flag & BAM_FDUP) && (mb->useDupes <= 1))
	return TRUE;
    else if ((core->flag & BAM_FDUP) && (mb->useDupes > 1))
    {
	uint8_t *zd_tag = bam_aux_get(bam, "ZD");
	if (zd_tag)
	{
	    char *s = bam_aux2Z(zd_tag);
	    char *num = chopPrefix(s);
	    unsigned n = sqlUnsigned(num);
	    if (n > mb->useDupes)
		return TRUE;
	}
    }
    if ((rg) && (mb->rgList) && ((!mb->rgListIsBlack && !hashLookup(mb->rgList, rg)) || (mb->rgListIsBlack && hashLookup(mb->rgList, rg))))
	return TRUE;
    return FALSE;
} 

static void metaBigAddBamFlagCounts(struct metaBig *mb, const bam1_core_t *core)
/* */
{
    int power;
    mb->bc.type_count[(int)core->flag]++;
    for (power = 0; power < 11; power++)
    {
	int decval = (int)pow(2,power);
	if (core->flag & decval)
	    mb->bc.bit_count[power]++;
    }
}

static int bamAddBed6(const bam1_t *bam, void *data)
/* bam_fetch() calls this on each bam alignment retrieved.  Translate each bam
 * into a bed. */
{
    struct metaBigBed6Helper *helper = (struct metaBigBed6Helper *)data;
    struct lm *lm = helper->lm;
    struct bed6 *bed;
    int start, end;
    struct metaBig *mb = helper->mb;
    const bam1_core_t *core = &bam->core;
    /* first find out if we're to skip this read because there's a strand choice */
    /* in the options */
    char strand = '+';
    /* don't consider if filtered */
    if (filterBam(mb, bam, core))
	return 0;
    if (bamIsRc(bam))
	strand = '-';
    /* check single or paired-end */
    if ((!(core->flag & BAM_FPAIRED)) || (core->flag & BAM_FMUNMAP))
    /* the read isn't paired or its mate is unmapped */ 
    {
	if ((mb->strand != 0) && (mb->strand != strand))
	    return 0;
	start = core->pos;
	/* do shifts, extensions */
	if (strand == '+')
	{
	    if (mb->length > 0)
		end = start + mb->length;
	    else if (core->isize > 0)
		end = start + core->isize;
	    else 
		end = bam_calend(core, bam1_cigar(bam));
	}
	else
	{
	    end = bam_calend(core, bam1_cigar(bam));
	    if (mb->length > 0)
		start = end - mb->length;
	    else if (core->isize > 0)
		start = end - core->isize;
	}
	/* filter out ones that go out-of-bounds after shifting/extending */
	if ((start >= end) || (start < 0) || (end > helper->chromSize))
	    return 0;
	/* ok fine, we're adding it to the bed. */
	lmAllocVar(lm, bed);
	metaBigAddBamFlagCounts(mb, core);
	if (mb->nameType == sequence)
	    bed->name = lmBamGetQuerySequence(bam, TRUE, lm);
	else if (mb->nameType == basicName)
	    bed->name = lmBamGetOriginalName(bam, lm);
	else if (mb->nameType == quality)
	bed->name = lmBamGetQuality(bam, TRUE, lm);
	else
	    bed->name = helper->dot;
	bed->chrom = helper->chrom;
	bed->chromStart = start;
	bed->chromEnd = end;
	bed->score = 1000;
	bed->strand[0] = strand;
	bed->strand[1] = '\0';
	slAddHead(&helper->bedList, bed);
    }
    else
    /* paired end */
    {
	if (core->flag & BAM_FPROPER_PAIR)
	{
	    /* only take one of the paired reads.  the other will be skipped */
	    if ((int)core->isize > 0)
	    {
		int mid;
		start = core->pos;
		end = start + core->isize;
		mid = core->pos + (core->isize / 2);
		if (mb->length > 0)
		{
		    start = mid - mb->length/2;
		    end = start + mb->length;
		}
		if ((start < 0) || (end > helper->chromSize))
		    return 0;
		lmAllocVar(lm, bed);
		metaBigAddBamFlagCounts(mb, core);
		bed->name = helper->dot;
		bed->chrom = helper->chrom;
		bed->chromStart = start;
		bed->chromEnd = end;
		bed->score = 1000;
		bed->strand[0] = '+';
		bed->strand[1] = '\0';
		slAddHead(&helper->bedList, bed);
	    }
	}
    }
    return 0;
}

static struct bed6 *bamBed6Fetch(struct metaBig *mb, char *chrom, unsigned start, unsigned end, struct lm *lm)
/* the main fetcher */
{
    struct metaBigBed6Helper helper;
    char posForBam[256];
    helper.lm = lm;
    helper.chrom = chrom;
    helper.dot = ".";
    helper.bedList = NULL;
    helper.mb = mb;
    helper.chromSize = hashIntVal(mb->chromSizeHash, chrom);
    safef(posForBam, sizeof(posForBam), "%s:%d-%d", chrom, start, end);
    bamFetchAlreadyOpen(mb->big.bam, mb->idx, mb->fileName, posForBam, bamAddBed6, &helper);
    slReverse(&helper.bedList);
    return helper.bedList;
}

static int bamAddPairbed(const bam1_t *bam, void *data)
/* bam_fetch() calls this on each bam alignment retrieved.  Translate each bam
 * into a bed. */
{
    struct metaBigPairbedHelper *helper = (struct metaBigPairbedHelper *)data;
    struct lm *lm = helper->lm;
    struct pairbed *pb;
    int start;
    char strand = '+';
    struct metaBig *mb = helper->mb;
    const bam1_core_t *core = &bam->core;
    if (filterBam(mb, bam, core))
	return 0;
    if (!((core->flag & BAM_FPAIRED) && (core->flag & BAM_FPROPER_PAIR)))
	return 0;
    if (!mb->useBothReads && (!(core->flag & BAM_FREAD1)))
	return 0;
    /* check whitelist/blacklist */
    if (bamIsRc(bam))
	strand = '-';
    start = core->pos;
    /* filter out ones that go out-of-bounds after shifting/extending */
    /* ok fine, we're adding it to the bed. */
    lmAllocVar(lm, pb);
    pb->name = lmBamGetOriginalName(bam, lm);
    pb->chrom = helper->chrom;
    pb->chromStart = start;
    pb->score = 1000;
    pb->strand[0] = strand;
    pb->strand[1] = '\0';
    pb->mChrom = helper->header->target_name[core->mtid];
    pb->mChromStart = core->mpos;
    slAddHead(&helper->pbList, pb);
    return 1;
}

static struct pairbed *bamPairbedFetch(struct metaBig *mb, char *chrom, unsigned start, unsigned end, struct lm *lm)
/* the main fetcher */
{
    struct metaBigPairbedHelper helper;
    char posForBam[256];
    helper.lm = lm;
    helper.chrom = chrom;
    helper.dot = ".";
    helper.pbList = NULL;
    helper.mb = mb;
    helper.header = mb->header;
    helper.chromSize = hashIntVal(mb->chromSizeHash, chrom);
    safef(posForBam, sizeof(posForBam), "%s:%d-%d", chrom, start, end);
    bamFetchAlreadyOpen(mb->big.bam, mb->idx, mb->fileName, posForBam, bamAddPairbed, &helper);
    slReverse(&helper.pbList);
    return helper.pbList;
}

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
	else
	    end += mb->shift;
    }
    else
    {
	if (mb->length > 0)
	    start -= (mb->shift + mb->length);
	else
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

static long bigBedCount(struct metaBig *mb, char *chrom, unsigned start, unsigned end)
/* main counter for bigBed.  I wonder if this can be done without fetching items */
{
    struct lm *lm = lmInit(0);
    struct bigBedInterval *ints = bigBedIntervalQuery(mb->big.bbi, chrom, (bits32)start, (bits32)end, 0, lm);
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

static void metaBigBamFlagCountsInit(struct metaBig *mb)
/* init the flags array */
{
    int i;
    for (i = 0; i < 11; i++)
	mb->bc.bit_count[i] = 0;
    for (i = 0; i < 2048; i++)
	mb->bc.type_count[i] = 0;
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
    if (mb->type == isaBam)
    {
	char *anotherFileName = NULL;
	mb->chromSizeHash = bamChromSizes(mb->fileName);
	mb->header = bamGetHeaderOnly(mb->fileName);
	mb->big.bam = bamOpen(mb->fileName, &anotherFileName);
	if (!sameString(mb->fileName, anotherFileName))
	    mb->fileName = anotherFileName;
	/* Also need to load the index since it's a bam */
	mb->idx = bam_index_load(mb->fileName);
    }
    else if (mb->type == isaBigBed)
    {
	mb->big.bbi = bigBedFileOpen(mb->fileName);
	mb->chromSizeHash = bbiChromSizes(mb->big.bbi);
    }
    else if (mb->type == isaBigWig)
    {
	mb->big.bbi = bigWigFileOpen(mb->fileName);
	mb->chromSizeHash = bbiChromSizes(mb->big.bbi);	
    }
    if (sectionsBed && sections)
    {
	struct bed *regions = bedLoadNAll(sectionsBed, 3);
	struct bed *subsets = subset_beds(sections, &regions, mb->chromSizeHash);
	mb->sections = subsets;
    }
    else if (sectionsBed)
    {
	mb->sections = bedLoadNAll(sectionsBed, 6);
    }
    else
	mb->sections = parseSectionString(sections, mb->chromSizeHash);
    metaBigBamFlagCountsInit(mb);
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
    if (mb->idx)
	bam_index_destroy(mb->idx);
    if (mb->type == isaBam)
	samclose(mb->big.bam);
    else if (mb->type == isaBigBed)
	bigBedFileClose(&mb->big.bbi);
    if (mb->header)
	bam_header_destroy(mb->header);
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
    if (mb->type == isaBam)
	return bamBed6Fetch(mb, chrom, start, end, lm);
    return NULL;
}

struct pairbed *metaBigPairbedFetch(struct metaBig *mb, char *chrom, unsigned start, unsigned end, struct lm *lm)
/* fetcher for pairbed for Hi-C */
{
    if (mb->type == isaBigBed)
	errAbort("bigBed not supported yet");
    if (mb->type == isaBam)
	return bamPairbedFetch(mb, chrom, start, end, lm);
    return NULL;
}

long metaBigCount(struct metaBig *mb, char *chrom, unsigned start, unsigned end)
/* the main counter */
{
    if (mb->type == isaBigBed)
	return bigBedCount(mb, chrom, start, end);
    if (mb->type == isaBam)
	return bamCount(mb, chrom, start, end);
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

static void binary_format(char *s, int num)
/* used but metaBigPrintFlagCounts to convert a number into a specific-length */
/* binary number. */
{
    int ix;
    for (ix = 0; ix < 11; ix++)
    {
	if ((int)pow(2,ix) & num)
	    s[10-ix] = '1';
	else 
	    s[10-ix] = '0';
    }
}

void metaBigPrintFlagCounts(struct metaBig *mb, char *file, boolean clear)
/* Print out the counts of the flags encountered with the last fetch function. */
{
    int i;
    char bin_num[12];
    FILE *output = NULL;
    if (file)
	output = mustOpen(file, "w");
    else 
	output = stdout;
    fprintf(output, "# flags:\n");
    fprintf(output, "    #  1    1 template having multiple segments in sequencing                           %'15ld\n", mb->bc.bit_count[0]);
    fprintf(output, "    #  2    2 each segment properly aligned according to the aligner                    %'15ld\n", mb->bc.bit_count[1]);
    fprintf(output, "    #  3    4 segment unmapped                                                          %'15ld\n", mb->bc.bit_count[2]);
    fprintf(output, "    #  4    8 next segment in the template unmapped                                     %'15ld\n", mb->bc.bit_count[3]);
    fprintf(output, "    #  5   16 SEQ being reverse complemented                                            %'15ld\n", mb->bc.bit_count[4]);
    fprintf(output, "    #  6   32 SEQ of the next segment in the template being reversed                    %'15ld\n", mb->bc.bit_count[5]);
    fprintf(output, "    #  7   64 the first segment in the template                                         %'15ld\n", mb->bc.bit_count[6]);
    fprintf(output, "    #  8  128 the last segment in the template                                          %'15ld\n", mb->bc.bit_count[7]);
    fprintf(output, "    #  9  256 secondary alignment                                                       %'15ld\n", mb->bc.bit_count[8]);
    fprintf(output, "    # 10  512 not passing quality controls                                              %'15ld\n", mb->bc.bit_count[9]);
    fprintf(output, "    # 11 1024 PCR or optical duplicate                                                  %'15ld\n", mb->bc.bit_count[10]);
    fprintf(output, "\n");
    bin_num[11] = '\0';
    if (clear)
	for (i = 0; i < 11; i++)
	    mb->bc.bit_count[i] = 0;
    for (i = 0; i < 2048; i++)
    {
	binary_format(bin_num, i);
	if (mb->bc.type_count[i] > 0)
    fprintf(output, " # %4d   %s         %'15ld\n", i, bin_num, mb->bc.type_count[i]);
	if (clear)
	    mb->bc.type_count[i] = 0;
    }
    if (file)
	carefulClose(&output);
}
