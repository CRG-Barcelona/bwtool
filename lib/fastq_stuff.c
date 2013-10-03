#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "common.h"
#include "hash.h"
#include "linefile.h"
#include "localmem.h"
#include "sqlNum.h"
#include "fastq_stuff.h"

boolean read_fastq_auto(struct fastq_auto *fq, struct lineFile *lf, boolean just_seq_qual)
/* fill in fastq struct from open lineFile.  return FALSE if at EOF */
/* set just_seq_qual=TRUE to skip loading everything except the sequence */
/* and quality information. */
{
    char *line;
    int len = 0;
    boolean neof = lineFileNext(lf, &line, &len);
    if (neof)
    {
	int i;
	int qual_size;
        /* should be header */
	if ((len <= 0) || (line[0] != '@'))
	    errAbort("Expecting header. Problem on line %d\n", lf->lineIx);
	if (!just_seq_qual)
	{
	    char *words[7];
	    int numWords;
	    numWords = chopByChar(line, ':', words, 6);
	    strcpy(fq->machine, words[0] + 1);
	    fq->flow_cell = sqlSigned(words[1]);
	    fq->tile = sqlSigned(words[2]);
	    fq->tile_x = sqlSigned(words[3]);
	    words[5] = chopPrefixAt(words[4], '#');
	    words[6] = chopPrefixAt(words[5], '/');
	    fq->tile_y = sqlSigned(words[4]);
	    fq->multiplex_index = sqlSigned(words[5]);
	    fq->pair_num = sqlSigned(words[6]);
	}
	/* read the sequence */
	fq->seq[0] = '\0';
	while ((neof = lineFileNext(lf, &line, &len)) && (len > 0) && (line[0] != '+'))
	    strcat(fq->seq, line);
	if (!neof)
	    errAbort("incomplete fastq file.  early EOF");
	fq->seq_size = strlen(fq->seq);
        /* at the point of the quality header.  who cares, read the quality */
	fq->qual[0] = '\0';
	while ((neof = lineFileNext(lf, &line, &len)) && (len > 0) && (line[0] != '@'))
	    strcat(fq->qual, line);
	if ((len > 0) && (line[0] == '@'))
	    lineFileReuse(lf);
	qual_size = strlen(fq->qual);
	if (qual_size != fq->seq_size)
	    errAbort("something wrong line %d.  sequence size (%d) should match quality size (%d)\n", 
		     lf->lineIx, fq->seq_size, qual_size);
	/* convert Illumina 1.3+ quals to Sanger */
	for (i = 0; i < qual_size; i++)
	    fq->qual[i] -= 64;
    }
    else 
	return FALSE;
    return TRUE;
}

void fprint_fastq_auto_normal(FILE *f, struct fastq_auto fq, int wrap_len)
/* print all the info from the fastq, essentially the same as it's read */
/* wrap_len specifies how long to make each line when printing the */
/* sequence and quality */
{
    fprintf(f, "@%s:%d:%d:%d:%d#%d/%d length=%d\n", fq.machine, fq.flow_cell, fq.tile, fq.tile_x, fq.tile_y, fq.multiplex_index, fq.pair_num, fq.seq_size);
    int i = 0;
    if (wrap_len > 0)
    {
	for (i = 0; i < fq.seq_size; i++)
	{
	    fprintf(f, "%c", fq.seq[i]);
	    if ((i+1 < fq.seq_size) && ((i+1) % wrap_len == 0))
		fprintf(f, "\n");
	}
    }
    fprintf(f, "\n+\n");
    if (wrap_len > 0)
    {
	for (i = 0; i < fq.seq_size; i++)
	{
	    fprintf(f, "%c", fq.qual[i] + 33);
	    if ((i+1 < fq.seq_size) && ((i+1) % wrap_len == 0))
		fprintf(f, "\n");
	}
    }
    fprintf(f, "\n");
}

void fastq_stats_init(struct fastq_stats *fqs)
/* init it all to zero even though it probably already is. */ 
{
    int i, j;
    fqs->total_count = 0;
    fqs->kept_count = 0;
    fqs->skipped_count = 0;
    fqs->longest_read = 0;
    fqs->all_b = 0;
    for (i = 0; i < FASTQ_MAX_LEN; i++)
    {
	fqs->original_lengths[i] = 0;
	fqs->trimmed_lengths[i] = 0;
	for (j = 0; j <= FASTQ_SANGER_MAX_QUAL; j++)
	{
	    fqs->before_quals[i][j] = 0;
	    fqs->after_quals[i][j] = 0;
	}
    }
}

void fastq_stats_write(char *filename, struct fastq_stats *fqs)
/* write it to a file */
{
    FILE *f = mustOpen(filename, "w");
    int i, j;
    if ((fqs->total_count > 0) && (fqs->longest_read > 0)) 
    {
	fprintf(f, "%lu\n%lu\n%lu\n%lu\n%u\n", fqs->total_count, fqs->all_b, fqs->skipped_count, fqs->kept_count, fqs->longest_read);
	for (i = 0; i < fqs->longest_read-1; i++)
	    fprintf(f, "%lu\t", fqs->original_lengths[i]);
	fprintf(f, "%lu\n", fqs->original_lengths[i]);
	for (i = 0; i < fqs->longest_read-1; i++)
	    fprintf(f, "%lu\t", fqs->trimmed_lengths[i]);
	fprintf(f, "%lu\n", fqs->trimmed_lengths[i]);
	for (j = 0; j <= FASTQ_SANGER_MAX_QUAL; j++)
	{
	    for (i = 0; i < fqs->longest_read-1; i++)
		fprintf(f, "%lu\t", fqs->before_quals[i][j]);
	    fprintf(f, "%lu\n", fqs->before_quals[i][j]);
	}
	for (j = 0; j <= FASTQ_SANGER_MAX_QUAL; j++)
	{
	    for (i = 0; i < fqs->longest_read-1; i++)
		fprintf(f, "%lu\t", fqs->after_quals[i][j]);
	    fprintf(f, "%lu\n", fqs->after_quals[i][j]);
	}
    }
    carefulClose(&f);	
}

void fastq_stats_read(char *filename, struct fastq_stats *fqs)
/* read from a file */
{
    struct lineFile *lf = lineFileOpen(filename, TRUE);
    char *line;
    char *words[FASTQ_MAX_LEN];
    int numWords = 0;
    int i, j;
    if (!lineFileNext(lf, &line, NULL))
	errAbort("bad file, total count section");
    fqs->total_count = sqlUnsignedLong(line);
    if (!lineFileNext(lf, &line, NULL))
	errAbort("bad file, skipped count section");
    fqs->all_b = sqlUnsignedLong(line);
    if (!lineFileNext(lf, &line, NULL))
	errAbort("bad file, skipped count section");
    fqs->skipped_count = sqlUnsignedLong(line);
    if (!lineFileNext(lf, &line, NULL))
	errAbort("bad file, skipped count section");
    fqs->kept_count = sqlUnsignedLong(line);
    if (!lineFileNext(lf, &line, NULL))
	errAbort("bad file, longest read section");
    fqs->longest_read = sqlUnsigned(line);
    numWords = lineFileChopTab(lf, words);
    if (numWords != fqs->longest_read)
	errAbort("bad file: original lengths section shoud have %u cols, has %d", fqs->longest_read, numWords);
    for (i = 0; i < numWords; i++)
	fqs->original_lengths[i] = sqlUnsignedLong(words[i]);
    numWords = lineFileChopTab(lf, words);
    if (numWords != fqs->longest_read)
	errAbort("bad file: trimmed lengths section shoud have %u cols, has %d", fqs->longest_read, numWords);
    for (i = 0; i < numWords; i++)
	fqs->trimmed_lengths[i] = sqlUnsignedLong(words[i]);
    for (i = 0; i <= FASTQ_SANGER_MAX_QUAL; i++)
    {
	numWords = lineFileChopTab(lf, words);
	if (numWords != fqs->longest_read)
	    errAbort("bad file:  before_quals section line %d shoud have %u cols, has %d", i+1, fqs->longest_read, numWords);
	for (j = 0; j < numWords; j++)
	    fqs->before_quals[j][i] = sqlUnsignedLong(words[j]);
    }
    for (i = 0; i <= FASTQ_SANGER_MAX_QUAL; i++)
    {
	numWords = lineFileChopTab(lf, words);
	if (numWords != fqs->longest_read)
	    errAbort("bad file:  after_quals section line %d shoud have %u cols, has %d", i+1, fqs->longest_read, numWords);
	for (j = 0; j < numWords; j++)
	    fqs->after_quals[j][i] = sqlUnsignedLong(words[j]);
    }
    lineFileClose(&lf);
}

void fastq_stats_add(struct fastq_stats *sum_stats, struct fastq_stats *fqs)
/* add all the numbers from the fqs one to the sum_stats */
{
    int i, j;
    sum_stats->total_count += fqs->total_count;
    sum_stats->skipped_count += fqs->skipped_count;
    sum_stats->all_b += fqs->all_b;
    sum_stats->kept_count += fqs->kept_count;
    if (fqs->longest_read > sum_stats->longest_read)
	sum_stats->longest_read = fqs->longest_read;
    for (i = 0; i < sum_stats->longest_read; i++)
    {
	sum_stats->original_lengths[i] += fqs->original_lengths[i];
	sum_stats->trimmed_lengths[i] += fqs->trimmed_lengths[i];
	for (j = 0; j <= FASTQ_SANGER_MAX_QUAL; j++)
	{
	    sum_stats->before_quals[i][j] += fqs->before_quals[i][j];
	    sum_stats->after_quals[i][j] += fqs->after_quals[i][j];
	}
    }
}
