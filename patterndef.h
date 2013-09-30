/*
 * patterndef.h
 *
 *  Created on: Apr 16, 2013
 *      Author: Houjun
 */
#include "uthash.h"
#ifndef PATTERNDEF_H
#define PATTERNDEF_H

#define CONTIGUOUS 1
#define SEQUENTIAL 2
#define KD_STRIDED 3

#define T_ADIO_READ 1
#define T_ADIO_WRITE 2
#define T_ADIO_READSTRIDED 3
#define T_ADIO_WRITESTRIDED 3

#define SEQ_LIST_THRESHOLD 3
#define PATTERN_SIZE_THRESHOLD 3
#define PATTERN_K_SIZE_MAX 10
#define MAX_LINE_LENGTH 256

typedef struct pattern {
    //char filepath[128];
    char operation;
    char patternType;
    int k;
    int mpiRank;
    uint64_t startPos;
    uint64_t  endPos;
    uint64_t  reqSize;
    uint64_t  strideSize[PATTERN_K_SIZE_MAX];
    uint64_t  recordNum[PATTERN_K_SIZE_MAX];
    double startTime;
    double endTime;
    uint64_t  reqOffesets[1024];
    struct pattern *prev;
    struct pattern *next;
} AccessPattern;

typedef struct tracelist {
    char op;
    int mpirank;
    uint64_t  offset;
    uint64_t  size;
    double startTime;
    double endTime;
    struct tracelist *prev;
    struct tracelist *next;
} TraceList;

// currently only use operation offset as key

typedef struct ut_lookup{
	uint64_t  key;	// key is offset
	AccessPattern *pattern;
    UT_hash_handle hh;
} UT_lookup;


typedef struct blocknode {
	int blocknum;		// key ~up to 2TB file with blksize 512B or more if size is flexible
	int freq;			// alignment? or use first few bits for size use & to mask
    UT_hash_handle hh;
} BlockNode;


typedef struct freq_block{
	int freq_blocknum;
	BlockNode *next;
    UT_hash_handle hh;
} Block_Hash;

// read from RADAR trace file
int read_radar(char* filename, TraceList** write_list, TraceList** read_list, AccessPattern** access_pattern
		, uint64_t * totalbytes);

// detect contiguous pattern from same process
int pattern_contig(TraceList** tracelist, AccessPattern** pattern_head, int mpirank);

// detect 1-D strided pattern from same process
int pattern_fixed_stride(TraceList** tracelist, AccessPattern** pattern_head, int mpirank);

// detect frequent access pattern
int freq_analysis(int totalcnt);

// lookup in hash table
int trace_lookup(TraceList **tracelist, TraceList *new_record, int mpirank);

// check if one trace record could be fitted into a known contiguous pattern
int contig_check(TraceList* trace, AccessPattern* pattern);

// check if one trace record could be fitted into a known 1-D strided pattern
int stride_check(TraceList* trace, AccessPattern* pattern);

// maps offset and lenth to block number, then sort
int offset_to_block(TraceList **tracelist);

// feed one trace record to each of the known patterns
int trace_feed(TraceList** tracelist, TraceList* new_record, AccessPattern** pattern_head, int mpirank);

// decide when to perform pattern detection algorithm to insert found pattern to known pattern list
int trace_analysis(TraceList** list, TraceList* tmp, AccessPattern** access_pattern, int* count, int mpirank);

// check if two pattern are same for patternType, strideSize, mpiRank and k value
int check_pattern_same(AccessPattern* p1, AccessPattern* p2);

// merge all (k-1)-D strided patterns to k-D pattern with given k
int merge_kd(AccessPattern** kdpattern_head, int k);

// allocate and add a trace to tracelist(a double linked list)
TraceList* addtmp(uint64_t  filepos, uint64_t  size, int op, double startTime, double endTime, int mpirank);

#endif
