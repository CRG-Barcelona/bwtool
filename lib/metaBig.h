/* metaBig - code to handle bigs the same way. particularly positional */
/* aspects of bigBed/bam.  tabix could be worked into here in the future */

#ifndef METABIG_H
#define METABIG_H

#ifndef HASH_H
#include "hash.h"
#endif

#ifndef BBIFILE_H
#include "bbiFile.h"
#endif

#ifndef BIGBED_H
#include "bigBed.h"
#endif

#ifdef USE_BAM
#include "bamFile.h"
#include "sam.h"
#include "bam.h"
#endif

enum metaBigFileType
/* basic enum to say what kind of file it is */
{
    isNotBig = 0,
    isaBigBed = 1,
    isaBigWig = 2,
    isaBam = 3,
    isaTabix = 4,
};

enum metaBigNameType
/* used for retrieval item-naming */
{
    justADot = 0,   /* (default) name with a '.' because often the names aren't needed */
    sequence = 1,   /* use sequence (if available) */
    basicName = 2,  /* use query name if BAM, or bed name if bed */
    quality = 3,    /* use sanger quality string (if available) */
};

struct bamFlagCounts
/* keep track of types of bams encountered */
{
    long type_count[2048];
    long bit_count[11];
};

struct metaBig
/* hold all the info we might ever need if it's a big */
{
    struct metaBig *next;          /* keep with sl-style structs */
    enum metaBigFileType type;     /* important thing to know what kinda file it is */
    union
    {                              /* should only be one of these things */
#ifdef USE_BAM
	samfile_t *bam;       
#endif
	struct bbiFile *bbi;
    } big;
    boolean isRemote;              /* is this a URL or not? */
    char *originalFileName;        /* full filename including URL or sectioning */
    char *fileName;                /* just the filename without sectioning */
    char *baseFileName;            /* the filename without remote dir or sectioning */
    char *remoteSiteAndDir;        /* the filename remote site+dir only.  NULL if remote==FALSE */
    struct bed *sections;          /* bed of sections from filename or an supplemented file */
    struct hash *chromSizeHash;    /* hash of chrom sizes taken from metadata in the big */
    /**** fetching options ****/
    int length;                    /* make all items n bases long (overrides extend) */
    int shift;                     /* shift items fetched by <shift> bases (can be negative) */
    char strand;                   /* restrict fetching to a particular strand, NULL if not */
    /* bam-only stuff */
    int mapQ;                      /* minimum mapping quality */
    long long numReads;                  /* number of reads total in the metaBig.  currently only set if using bigBed */
    boolean pe;                    /* only use paired-end reads */
    boolean se;                    /* only use single-end reads */
    boolean useMapQ;               /* use the MAPQ setting... otherwise rely on flags in the bam for QC fail */
    boolean unpaired;              /* return reads without pairing them */
    boolean includeB;              /* include B-reads */
    boolean includeBadRegions;     /* include bad regions */
    boolean includeMissingRestrSite;
    int useDupes;                  /* */
    boolean useBothReads;          /* for paired-bed retrieval */
    enum metaBigNameType nameType; /* the type of naming given to each item fetched */
    struct hash *rgList;           /* just a list of read-groups to restrict to (a black or whitelist)  */
    boolean rgListIsBlack;         /* by default the list is a whitelist, otherwise it's a blacklist if */
                                   /* this is TRUE. */
#ifdef USE_BAM
    bam_index_t *idx;              /* NULL if not a BAM file */ 
    struct bamFlagCounts bc;       /* counts of the flags of the used reads */ 
    bam_header_t *header;          /* header info */
#endif
};

struct pairbed
{
    struct pairbed *next;
    char *chrom;
    unsigned chromStart;
    char *name;
    int score;
    char strand[2];
    char *mChrom;
    unsigned mChromStart;
};

struct bed6
/* just the first 6 fields of a bed. for high-mem situations, it's a little */
/* nicer not to have the whole enchilada of bed 15 */
{
    struct bed6 *next;  /* Next in singly linked list. */
    char *chrom;	/* Human chromosome or FPC contig */
    int chromStart;	/* Start position in chromosome */
    int chromEnd;	/* End position in chromosome */
    char *name;	/* Name of item */
    int score; /* Score - 0-1000 */
    char strand[2];  /* + or -.  */    
};

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

#ifdef USE_BAM
#include "metaBigBam.h"
#endif

struct bed6 *readBed6(char *file);
/* read from a file */

struct bed6 *readBed6Soft(char *file);
/* read from a file.  If it's missing fields, fill the bed */

void bed6Free(struct bed6 **pEl);
/* Free a single dynamically allocated bed such as created
 * with bedLoad(). */

void bed6FreeList(struct bed6 **pList);
/* Free a list of dynamically allocated bed's */

enum metaBigFileType isBigWigOrBed(char *filename);
/* Peak at a file to see if it's bigWig */

enum metaBigFileType sniffBigFile(char *filename);
/* try to figure out what type of big file it is. */

struct metaBig *metaBigOpen(char *fileOrUrlwSections, char *sectionsBed);
/* load a file or URL with or without sectioning */
/* if it's a bam, load the index. */

void metaBigSetPositionalOptions(struct metaBig *mb, int length, int shift, char strand);
/* set retrieval options for the metaBig. */

void metaBigSetNameOption(struct metaBig *mb, enum metaBigNameType);
/* set the naming option for retrieval */

void metaBigSetRgList(struct metaBig *mb, char *commaSep, boolean isBlackList);
/* sets the whitelist (or blacklist) */

void metaBigToggleBadRegions(struct metaBig *mb, boolean useBad);
/* setting this to TRUE permits the retrieval of reads marked as from "bad regions" */
/* from the metaBig (only relevant with BAM files using the ZR tag) */

void metaBigToggleBRead(struct metaBig *mb, boolean useB);
/* setting this to TRUE permits the retrieval of reads marked as B-reads from the metaBig */
/* (only relevant with BAM files using the ZL tag) */

void metaBigClose(struct metaBig **pMb);
/* close the file and free up everything. */

struct bed6 *metaBigBed6Fetch(struct metaBig *mb, char *chrom, unsigned start, unsigned end, struct lm *lm);
/* the main fetcher */

struct pairbed *metaBigPairbedFetch(struct metaBig *mb, char *chrom, unsigned start, unsigned end, struct lm *lm);
/* fetcher for pairbed for Hi-C */

long metaBigCount(struct metaBig *mb, char *chrom, unsigned start, unsigned end);
/* the main counter */

long metaBigNumItems(struct metaBig *mb, boolean verbose);
/* return the total number of items in a bigBed or BAM */
/* used on a bigWig will return 0 */
/* unfortunately this is a loop through the entire file basically. */
/* nicer would be something that just glances at the index, but doing that */
/* might count items that would be filtered out upon fetching. */

void metaBigPrintFlagCounts(struct metaBig *mb, char *file, boolean clear);
/* Print out the counts of the flags encountered with the last fetch function. */

#endif /* METABIG_H */
