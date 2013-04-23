#ifndef FASTQ_STUFF_H
#define FASTQ_STUFF_H

#define FASTQ_MAX_LEN 2048
#define FASTQ_SANGER_MAX_QUAL 40

struct fastq_auto
/* fastq using all non-heap mem */
{
    char seq[FASTQ_MAX_LEN];
    char qual[FASTQ_MAX_LEN];
    int seq_size;
    char machine[128];
    int flow_cell;
    int tile;
    int tile_x;
    int tile_y;
    int multiplex_index;
    int pair_num;
};

boolean read_fastq_auto(struct fastq_auto *fq, struct lineFile *lf, boolean just_seq_qual);
/* fill in fastq struct from open lineFile.  return FALSE if at EOF */
/* set just_seq_qual=TRUE to skip loading everything except the sequence */
/* and quality information. */

void fprint_fastq_auto_normal(FILE *f, struct fastq_auto fq, int wrap_len);
/* print all the info from the fastq, essentially the same as it's read */
/* wrap_len specifies how long to make each line when printing the */
/* sequence and quality */

struct fastq_stats
{
    long unsigned before_quals[FASTQ_MAX_LEN][FASTQ_SANGER_MAX_QUAL+1];
    long unsigned after_quals[FASTQ_MAX_LEN][FASTQ_SANGER_MAX_QUAL+1];
    long unsigned total_count;
    long unsigned all_b;
    long unsigned skipped_count;
    long unsigned kept_count;
    long unsigned original_lengths[FASTQ_MAX_LEN];
    long unsigned trimmed_lengths[FASTQ_MAX_LEN];
    unsigned longest_read;
};

void fastq_stats_init(struct fastq_stats *fqs);
/* init it all to zero even though it probably already is. */ 

void fastq_stats_write(char *filename, struct fastq_stats *fqs);
/* write it to a file */

void fastq_stats_read(char *filename, struct fastq_stats *fqs);
/* read from a file */

void fastq_stats_add(struct fastq_stats *sum_stats, struct fastq_stats *fqs);
/* add all the numbers from the fqs one to the sum_stats */

#endif
