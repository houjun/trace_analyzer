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

#define LIST_SIZE_THRESHOLD 3
#define PATTERN_SIZE_THRESHOLD 3
#define PATTERN_K_SIZE_MAX 10
#define MAX_LINE_LENGTH 256

typedef struct pattern {
    //char filepath[128];
    char operation;
    char patternType;
    int k;
    int mpiRank;
    int startPos;
    int endPos;
    int reqSize;
    int strideSize[PATTERN_K_SIZE_MAX];
    int recordNum[PATTERN_K_SIZE_MAX];
    double startTime;
    int reqOffesets[1024];
    struct pattern *prev;
    struct pattern *next;
} AccessPattern;


typedef struct tracelist {
    char op;
    int mpirank;
    int offset;
    int size;
    double opTime;
    struct tracelist *prev;
    struct tracelist *next;
} TraceList;

// currently only use operation, rank, offset as key
// request size not used because of sequential pattern
typedef struct lookup_key{
	char op;
	int mpirank;
	int off;
} Lookup_key;

typedef struct {
	Lookup_key key;
	AccessPattern *pattern;
    UT_hash_handle hh;
}UT_lookup;


// read from RADAR trace file
int read_radar(char* filename, TraceList** write_list, TraceList** read_list, AccessPattern** access_pattern
		, int* totalbytes, UT_lookup **lookup);

// detect contiguous pattern from same process
int pattern_contig(TraceList** tracelist, AccessPattern** pattern_head, int mpirank, UT_lookup **lookup);

// detect 1-D strided pattern from same process
int pattern_fixed_stride(TraceList** tracelist, AccessPattern** pattern_head, int mpirank, UT_lookup **lookup);

// check if one trace record could be fitted into a known contiguous pattern
int contig_check(TraceList* trace, AccessPattern* pattern);

// check if one trace record could be fitted into a known 1-D strided pattern
int stride_check(TraceList* trace, AccessPattern* pattern);

// feed one trace record to each of the known patterns
int trace_feed(TraceList** tracelist, TraceList* new_record, AccessPattern** pattern_head, int mpirank);

// decide when to perform pattern detection algorithm to insert found pattern to known pattern list
int trace_analyzer(TraceList** list, TraceList* tmp, AccessPattern** access_pattern, int* count, int mpirank
		, UT_lookup **lookup);

// check if two pattern are same for patternType, strideSize, mpiRank and k value
int check_pattern_same(AccessPattern* p1, AccessPattern* p2);

// merge all (k-1)-D strided patterns to k-D pattern with given k
int merge_kd(AccessPattern** kdpattern_head, int k);

// allocate and add a trace to tracelist(a double linked list)
TraceList* addtmp(int filepos, int size, int op, double opTime, int mpirank);

#endif
