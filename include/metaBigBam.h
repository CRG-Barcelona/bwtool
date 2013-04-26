/* DON'T include this header directly, rather metaBig.h */
/* This one is only meant to have some bam-related code */

#ifndef METABIGBAM_H
#define METABIGBAM_H

boolean isBamWithIndex(char *file);
/* attempt to open a BAM file and its index. if everything is ok, */
/* return TRUE before closing. */

bam_header_t *bamGetHeaderOnly(char *file);
/* just get the header, close the BAM */

struct hash *bamChromSizes(char *bamfile);
/* return a list of chromosome sizes from a bam file */

void metaBigBamFlagCountsInit(struct metaBig *mb);
/* init the flags array */

long bamCount(struct metaBig *mb, char *chrom, unsigned start, unsigned end);
/* the main fetcher of counts of bams */

struct bed6 *bamBed6Fetch(struct metaBig *mb, char *chrom, unsigned start, unsigned end, struct lm *lm);
/* the main fetcher of bams */

struct pairbed *bamPairbedFetch(struct metaBig *mb, char *chrom, unsigned start, unsigned end, struct lm *lm);
/* the main fetcher of bams as pairs */

void metaBigPrintFlagCounts(struct metaBig *mb, char *file, boolean clear);
/* Print out the counts of the flags encountered with the last fetch function. */

#endif  /* METABIGBAM_H */
