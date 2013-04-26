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

#ifdef USE_BAM

/* functions for bam-related things, which are separate here in case samtools isn't linked */

/**** file-sniffing code ****/

typedef struct {  
    int beg, end;  
    samfile_t *in;  
} tmpstruct_t;

boolean isBamWithIndex(char *file)
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

bam_header_t *bamGetHeaderOnly(char *file)
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

struct hash *bamChromSizes(char *bamfile)
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

void metaBigBamFlagCountsInit(struct metaBig *mb)
/* init the flags array */
{
    int i;
    for (i = 0; i < 11; i++)
	mb->bc.bit_count[i] = 0;
    for (i = 0; i < 2048; i++)
	mb->bc.type_count[i] = 0;
}

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

long bamCount(struct metaBig *mb, char *chrom, unsigned start, unsigned end)
/* the main fetcher of counts of bams */
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

struct bed6 *bamBed6Fetch(struct metaBig *mb, char *chrom, unsigned start, unsigned end, struct lm *lm)
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

struct pairbed *bamPairbedFetch(struct metaBig *mb, char *chrom, unsigned start, unsigned end, struct lm *lm)
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
#endif /* USE_BAM */
