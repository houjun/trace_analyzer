/*
 * patterndef.h
 *
 *  Created on: Apr 16, 2013
 *      Author: Houjun
 */
#ifndef PATTERNDEF_H
#define PATTERNDEF_H

#define CONTIGUOUS 1
#define REPEAT 2
#define KD_STRIDED 3
#define SEQUENTIAL 4

#define LIST_SIZE_THRESHOLD 3
#define PATTERN_SIZE_THRESHOLD 3
#define PATTERN_K_SIZE_MAX 10
#define MAX_LINE_LENGTH 256

/*******************************************
 * Supported trace operations (from trace.h)
 *******************************************/
#define T_ADIO_OPEN           1
#define T_ADIO_OPENCOLL       2
#define T_ADIO_CLOSE          3

#define T_ADIO_READCONTIG     4
#define T_ADIO_READSTRIDED    5
#define T_ADIO_READCOLL       6
#define T_ADIO_WRITECONTIG    7
#define T_ADIO_WRITESTRIDED   8
#define T_ADIO_WRITECOLL      9

#define T_ADIO_IREADCONTIG   10
#define T_ADIO_IWRITECONTIG  11
#define T_ADIO_IREADSTRIDED  12
#define T_ADIO_IWRITESTRIDED 13

#define T_ADIO_FCNTL         14
#define T_ADIO_FEATURE       15
#define T_ADIO_FLUSH         16
#define T_ADIO_RESIZE        17
#define T_ADIO_SEEKIND       18
#define T_ADIO_SETINFO       19
//above 19 are unused

#define T_ADIO_READ 1
#define T_ADIO_WRITE 2

typedef struct pattern {
    //char filepath[128];
    char operation;
    char patternType;
    int k;
    int mpiRank;
    int startPos;
    int endPos;
    int startStep;
    int endStep;
    int reqSize;
    int strideSize[PATTERN_K_SIZE_MAX];
    int recordNum[PATTERN_K_SIZE_MAX];	// # of records
    double startTime;
    int reqOffesets[1024];
    struct pattern *prev;
    struct pattern *next;
} AccessPattern;

typedef struct markovpattern {
    //char filepath[128];
    char operation;
    char patternType;
    int k;
    int mpiRank;
    int startPos[PATTERN_K_SIZE_MAX];
    int endPos[PATTERN_K_SIZE_MAX];
    int reqSize[PATTERN_K_SIZE_MAX];
    int strideSize[PATTERN_K_SIZE_MAX];
    int recordNum[PATTERN_K_SIZE_MAX];
    int startStep;
    int endStep;
    double startTime;
    struct markovpattern *prev;
    struct markovpattern *next;
} Markovpattern;
//Markov pattern not used for now

typedef struct tracelist {
    char op;
    int mpirank;
    int offset;
    int size;
    int opStep;
    double opTime;
    struct tracelist *prev;
    struct tracelist *next;
} TraceList;

// read from RADAR trace file
int read_radar(char* filename, TraceList** write_list, TraceList** read_list, AccessPattern** access_pattern, int* totalbytes);

// detect contiguous pattern from same process
int pattern_contig(TraceList** tracelist, AccessPattern** pattern_head, int mpirank);

// detect 1-D strided pattern from same process
int pattern_fixed_stride(TraceList** tracelist, AccessPattern** pattern_head, int mpirank);

// detect markov pattern from same process, not used for now
int pattern_markov(TraceList** tracelist, AccessPattern** pattern_head, int mpirank);

// check if one trace record could be fitted into a known contiguous pattern
int contig_check(TraceList* trace, AccessPattern* pattern);

// check if one trace record could be fitted into a known 1-D strided pattern
int stride_check(TraceList* trace, AccessPattern* pattern);

// check if one trace record could be fitted into a known repeat pattern
int repeat_check(TraceList* trace, AccessPattern* pattern);

// feed one trace record to each of the known patterns
int trace_feed(TraceList** tracelist, TraceList* new_record, AccessPattern** pattern_head, int mpirank);

// decide when to perform pattern detection algorithm to insert found pattern to known pattern list
int trace_analyzer(TraceList** list, TraceList* tmp, AccessPattern** access_pattern, int* count, int mpirank);

// check if two pattern are same for patternType, strideSize, mpiRank and k value
int check_pattern_same(AccessPattern* p1, AccessPattern* p2);

// merge markov pattern to compact format
int merge_markov(AccessPattern** pattern_head, Markovpattern** kdstride_head);

// merge all (k-1)-D strided patterns to k-D pattern with given k
int merge_kd(AccessPattern** kdpattern_head, int k);

// allocate and add a trace to tracelist(a double linked list)
TraceList* addtmp(int filepos, int size, int op, int opStep, double opTime, int mpirank);

#endif
