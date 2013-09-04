/*
 *pattern.c
 *
 * Created on: Apr 16, 2013
 *     Author: Houjun
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "patterndef.h"
#include "utlist.h"
#include "uthash.h"

extern int commsize;
extern int block_size;
extern int candidate_freq;
extern int result_freq;
extern int lookahead_window;
extern int freq_list_size;
extern int *count_nonseq;
extern int *block_index;
extern int block_index_size;
extern int *block_num;
extern int block_num_size;
extern Block_Hash *candidate;
extern Block_Hash *result;

extern TraceList *tmp_record;
extern UT_lookup lookup_tmp;
extern UT_lookup *found;
extern UT_lookup *lookup_replace;
extern UT_lookup *lookup;
int read_radar(char *filename, TraceList **write_list, TraceList **read_list, AccessPattern **access_pattern, int *totalbytesrw)
{

	char buf[MAX_LINE_LENGTH];
	FILE *fp = NULL;
	if( (fp = fopen(filename, "r+")) == NULL) {
	    printf("No such file\n");
	    return -1;
	}
	if (fp == NULL)	{
	    printf("Error Reading File\n");
	}

	int i, j, nrecords, mpairs;
	char operator[MAX_LINE_LENGTH];
	char filepath[MAX_LINE_LENGTH];
	char processtmp[MAX_LINE_LENGTH];
	int read_count = 0;
	int write_count = 0;
	int hierarchy_info;
	int proc_num;

	while(fgets(buf, MAX_LINE_LENGTH, fp)){
		if (*buf == '#') continue; // ignore comment line

		if(strcmp(buf,"HEADER\n")==0){
			//read next lines as header
			//<Full File path> <MPI file comm size>
			fgets(buf, MAX_LINE_LENGTH, fp);
			sscanf(buf, "%s %d %d",filepath, &commsize, totalbytesrw);
			proc_num = commsize;

			while(proc_num--){
				//PROCESS <MPI Rank> <n: number of trace entries>
				fgets(buf, MAX_LINE_LENGTH, fp);
				sscanf(buf, "%s %d %d",processtmp, &tmp_record->mpirank, &nrecords);
				for(i=0; i < nrecords; i++){

					//read next n records
					fgets(buf, MAX_LINE_LENGTH, fp);
					//<Op string> <Time delta> <m: # of accesses>
					sscanf(buf, "%d %s %lf %lf %d",&hierarchy_info, operator, &tmp_record->startTime, &tmp_record->endTime, &mpairs);

					for(j=0; j < mpairs; j++){
						fgets(buf, MAX_LINE_LENGTH, fp);
						sscanf(buf, "%d, %d", &tmp_record->offset, &tmp_record->size);

						//add trace record
						if(strstr(operator,"Read") != NULL){
							if(strstr(operator,"Strided") != NULL)
								tmp_record->op = T_ADIO_READSTRIDED;
							else
								tmp_record->op = T_ADIO_READ;

							trace_analysis(read_list, tmp_record, access_pattern, &read_count, tmp_record->mpirank);
							//just like below
						}
						else if(strstr(operator,"Write") != NULL){
							if(strstr(operator,"Strided") != NULL)
								tmp_record->op = T_ADIO_WRITESTRIDED;
							else
								tmp_record->op = T_ADIO_WRITE;

							trace_analysis(write_list, tmp_record, access_pattern, &write_count, tmp_record->mpirank);
						}

					}

				}

				//perform final analysis when reach end of each process
				trace_analysis(read_list, NULL, access_pattern, &read_count, tmp_record->mpirank);
				trace_analysis(write_list, NULL, access_pattern, &write_count, tmp_record->mpirank);


				read_count = 0;
				write_count = 0;

				TraceList *elt;
				TraceList *etmp;
				//clean up to analyze next process's trace records
				DL_FOREACH_SAFE(*read_list,elt,etmp) {
				      DL_DELETE(*read_list,elt);
				}
				DL_FOREACH_SAFE(*write_list,elt,etmp) {
					  DL_DELETE(*write_list,elt);
				}
				// clean up hash table
				UT_lookup *current_lookup, *tmp;
				HASH_ITER(hh, lookup, current_lookup, tmp) {
					HASH_DEL(lookup,current_lookup);
					free(current_lookup);
				}

			}

		}//Header

	}

	fclose(fp);
	return 0;
}

int trace_analysis(TraceList **list, TraceList *tmp, AccessPattern **access_pattern, int *count, int mpirank)
{
	int decreasecount;

	if(tmp != NULL && !trace_lookup(list, tmp, mpirank)){
		(*count)++;
		TraceList *add_tmp = addtmp(tmp->offset, tmp->size, tmp->op, tmp->startTime, tmp->endTime, tmp->mpirank);
		DL_APPEND(*list, add_tmp);
	}
	else
		return 0;

	// perform sequential analysis
	if(*count > SEQ_LIST_THRESHOLD || tmp == NULL){
		decreasecount = pattern_contig(list, access_pattern, mpirank);
		decreasecount = pattern_fixed_stride(list, access_pattern, mpirank);

		*count = *count - decreasecount;
	}

	/*
	// freq analysis
	if(*count > freq_list_size || tmp == NULL){

		// map to block
		int cnt = offset_to_block(list);
		freq_analysis(cnt);

		*count = 0;
	}
*/
	return 0;
}


int trace_lookup(TraceList **tracelist, TraceList *new_record, int mpirank)
{

	int findpattern = 0;

	lookup_tmp.key.op = new_record->op;
	lookup_tmp.key.off = new_record->offset;
	lookup_tmp.key.mpirank = new_record->mpirank;

	// look up
	HASH_FIND(hh, lookup, &lookup_tmp.key, sizeof(Lookup_key), found);

	// found corresponding pattern
	if(found != NULL){
		// update pattern and hash table
		if(found->pattern->patternType== CONTIGUOUS || found->pattern->patternType == SEQUENTIAL){

			lookup_replace = (UT_lookup*)malloc( sizeof(UT_lookup) );
			lookup_replace->pattern = found->pattern;
			lookup_replace->pattern->recordNum[0] ++;
			lookup_replace->pattern->endPos = new_record->offset;
			lookup_replace->key.op= found->key.op;
			lookup_replace->key.mpirank = found->key.mpirank;
			lookup_replace->pattern->reqSize = new_record->size;

			if(found->pattern->patternType == SEQUENTIAL){
				lookup_replace->pattern->reqOffesets[lookup_replace->pattern->recordNum[0]-1] = new_record->size;
			}

			lookup_replace->key.off = new_record->offset + new_record->size;
			HASH_DEL(lookup, found);
		    HASH_ADD(hh, lookup, key, sizeof(Lookup_key), lookup_replace);

			free(found);
			findpattern++;
		}
		else if(found->pattern->patternType == KD_STRIDED){

			lookup_replace = (UT_lookup*)malloc( sizeof(UT_lookup) );
			lookup_replace->pattern = found->pattern;
			lookup_replace->pattern->recordNum[0] ++;
			lookup_replace->pattern->endPos = new_record->offset;
			lookup_replace->key.op= found->key.op;
			lookup_replace->key.mpirank = found->key.mpirank;
			lookup_replace->key.off = new_record->offset + found->pattern->strideSize[0];

			//HASH_REPLACE(hh, *lookup, key, sizeof(Lookup_key), lookup_replace, found);
			HASH_DEL(lookup, found);
		    HASH_ADD(hh, lookup, key, sizeof(Lookup_key), lookup_replace);

			free(found);
			findpattern++;
		}
	}

	if(findpattern == 0){
		if(new_record==NULL){
			printf("Can't feed NULL!\n");
			return 0;
		}

	}

	return findpattern;
}

int sort_by_freq(BlockNode *a, BlockNode *b)
{
	 if (a->freq == b->freq) return 0;
	 return (a->freq < b->freq) ? -1 : 1;
}

int compareBlockNum(const void * a, const void * b)
{
	return ( *(int*)a - *(int*)b );
}

int offset_to_block(TraceList **tracelist)
{
	// block_index stores all unique block numbers
	// count_nonseq stores corresponding counts

	TraceList *list_i = *tracelist, *old;
	int i = 0, j = 0, totalcnt;

	while(list_i != NULL){
		// skip write op
		if(list_i->op % 2 == 0){
			list_i = list_i->next;
			continue;
		}

		for(j = 0; j < list_i->size / block_size; j++){
			if(i >= block_num_size)
			block_num = realloc(block_num, sizeof(int) * block_num_size * 2);
			block_num[i++] = list_i->offset / block_size;
		}
		old = list_i;
		list_i = list_i->next;
		DL_DELETE(*tracelist, old);

	}

	totalcnt = i;
	qsort(block_num, totalcnt, sizeof(int), compareBlockNum);

	j = 1;
	block_index[0] = block_num[0];
	for(i = 1; i < totalcnt; i++){
		if(j >= block_index_size){
			count_nonseq = realloc(count_nonseq, sizeof(int) * block_index_size * (totalcnt / i));
			block_index = realloc(block_index, sizeof(int) * block_index_size * (totalcnt / i));
		}

		if(block_num[i] == block_num[i - 1]){
			count_nonseq[j]++;
		}
		else{
			block_index[j++] = block_num[i];
		}

	}

	return totalcnt;
}

int freq_analysis(int totalcnt)
{
	int i, j;

	for(i = 0; i < totalcnt; i++){
		// first check if current(block_index[i]) exists in result or candidate hash table
		// if so update, else insert

		Block_Hash *tmp_candidate = NULL, *tmp_result = NULL;
		Block_Hash *tmp_switch = NULL;
		BlockNode *tmp_blk = NULL;

		HASH_FIND_INT(candidate, &block_index[i], tmp_candidate);
		HASH_FIND_INT(result, &block_index[i], tmp_result);

		if(tmp_candidate != NULL || tmp_result != NULL){
			if(tmp_result != NULL)
				tmp_switch = tmp_result;
			else
				tmp_switch = tmp_candidate;
		}
		else if(count_nonseq[i] > candidate_freq){
			// add at least to candidate set
			if(count_nonseq[i] > result_freq)
				tmp_switch = tmp_result;
			else
				tmp_switch = tmp_candidate;

			Block_Hash *s = malloc(sizeof(Block_Hash));
			s->freq_blocknum = block_index[i];
			s->next = NULL;

			HASH_ADD_INT(tmp_switch, freq_blocknum, s );

		}
		else{
			continue;
		}
		// next lookahead_window accesses
		// careful with boundary
		for(j = i + 1; j < i + lookahead_window && j < freq_list_size; j++){
			HASH_FIND_INT(tmp_switch->next, &block_index[j], tmp_blk);

			if(tmp_blk != NULL){
				// found previous one, just update
				tmp_blk->freq++;
				continue;
			}
			tmp_blk = malloc(sizeof(BlockNode));
			tmp_blk->blocknum = block_index[j];
			tmp_blk->freq = 1;

			HASH_ADD_INT(tmp_switch->next, blocknum, tmp_blk);

		}

	}//for
	return 0;
}


int pattern_repeat(TraceList **tracelist, AccessPattern **pattern_head, int mpirank, UT_lookup **lookup)
{
	TraceList *list_i = *tracelist;
	//TraceList *list_j = NULL;

	while(list_i != NULL && list_i->next != NULL){

	}

	return 0;
}

TraceList* addtmp(int filepos, int size, int op, double startTime, double endTime, int mpirank){

	TraceList *tmp ;
	if ( (tmp = (TraceList*)malloc(sizeof(TraceList))) == NULL)
		return NULL;
	tmp->mpirank = mpirank;
	tmp->offset = filepos;
	tmp->size = size;
	tmp->op = op;
	tmp->startTime = startTime;
	tmp->endTime = endTime;
	tmp->next = NULL;
	tmp->prev = NULL;
	return tmp;
}

int contig_check(TraceList *trace, AccessPattern *pattern)
{
	//TODO support for same offset but different size
	if(trace->mpirank == pattern->mpiRank && trace->offset == pattern->endPos + pattern->reqSize)
		return 1;
	else
		return 0;
}

int stride_check(TraceList *trace, AccessPattern *pattern)
{
	//TODO support for same offset but different size
	if(trace->mpirank == pattern->mpiRank && trace->offset == pattern->endPos  + pattern->strideSize[0]
	         && trace->size == pattern->reqSize)
		return 1;
	else
		return 0;
}

int pattern_contig(TraceList **tracelist, AccessPattern **pattern_head, int mpirank)
{
	TraceList *list_i = *tracelist;
	int isseq = 0;

	while(list_i != NULL && list_i->next != NULL){

		int contig_size = 1;
		isseq = 0;
		int req_arr_size = 0;

		//check each other records in the list to see if forms contig pattern
		TraceList *list_j = list_i->next;
		int tmppos = list_i->offset;
		int tmpsize = list_i->size;


		while(list_j != NULL){
			if(list_j->offset == tmppos + tmpsize){
				contig_size ++;
				tmppos = tmppos + tmpsize;
				if(list_j->size != tmpsize)
					isseq = 1;
			}
			list_j = list_j->next;


			if(contig_size >= PATTERN_SIZE_THRESHOLD){
				// add new pattern
				AccessPattern *contig_pattern ;
				if ( (contig_pattern = (AccessPattern*)malloc(sizeof(AccessPattern)) ) == NULL){
					printf("Failed to allocate space for Access Pattern!");
					return -1;
				}

				if(isseq == 0)
					contig_pattern->patternType = CONTIGUOUS;
				else
					contig_pattern->patternType = SEQUENTIAL;

				contig_pattern->k = 1;
				contig_pattern->reqSize = list_i->size;
				contig_pattern->startPos = list_i->offset;
				contig_pattern->endPos = list_i->offset;
				contig_pattern->operation = list_i->op;
				contig_pattern->strideSize[0] = 0;
				contig_pattern->recordNum[0] = 1;
				contig_pattern->mpiRank = mpirank;
				contig_pattern->startTime = list_i->startTime;
				contig_pattern->endTime = list_i->endTime;
				contig_pattern->reqOffesets[req_arr_size++] = list_i->size;

				DL_APPEND(*pattern_head, contig_pattern);

				// Delete the records corresponding to this pattern, feed the rest to this pattern
				TraceList *old = list_i;
				TraceList *oldnext = old->next;
				DL_DELETE(*tracelist, old);

				old = oldnext;
				oldnext = old->next;
				while(1){
					if(contig_check(old, contig_pattern)){
						// old in the pattern
						contig_pattern->endPos = old->offset;
						contig_pattern->reqSize = old->size;
						contig_pattern->recordNum[0] ++;
						contig_pattern->reqOffesets[req_arr_size++] = old->size;

						DL_DELETE(*tracelist, old);

					}
					old = oldnext;
					if(old == NULL)
						break;
					oldnext = old->next;
				}

				// add to look up hash table
				UT_lookup *tmp = (UT_lookup*)malloc( sizeof(UT_lookup) );
			    memset(tmp, 0, sizeof(UT_lookup));
				tmp->key.op = contig_pattern->operation;
				tmp->key.mpirank = contig_pattern->mpiRank;
				tmp->key.off = contig_pattern->endPos + contig_pattern->reqSize;
				tmp->pattern = contig_pattern;

			    HASH_ADD(hh, lookup, key, sizeof(Lookup_key), tmp);

				// if one contig pattern is found then return
				return contig_pattern->recordNum[0];
			}

		}

		list_i = list_i->next;

	}

	return 0;
}

int pattern_fixed_stride(TraceList **tracelist, AccessPattern **pattern_head, int mpirank)
{
	TraceList *list_i = *tracelist;
	TraceList *list_j = NULL;
	TraceList *list_k = NULL;

	while(list_i != NULL && list_i->next != NULL){

		list_j = list_i->next;
		int stride_size = 2;
		int tmppos = 0;
		int tmp_stride_size = 0;

		while(list_j != NULL){
			if(list_i->size != list_j->size){
				list_j = list_j->next;
				continue;
			}
			tmppos = list_j->offset;
			tmp_stride_size = list_j->offset - list_i->offset;
			// stride size should be greater than 0
			if(tmp_stride_size <= 0)
				break;

			list_k = list_j->next;

			//check each pair of records form a stride that the rest of records from the list could fit in
			while(list_k != NULL ){
				if(list_k->size != list_j->size){
					list_k = list_k->next;
					continue;
				}
				if(list_k->offset == tmppos + tmp_stride_size){
					stride_size ++;
					tmppos = tmppos + tmp_stride_size;
				}
				list_k = list_k->next;
			}

			if(stride_size >= PATTERN_SIZE_THRESHOLD){
				// add new pattern
				AccessPattern *stride_pattern ;
				if ( (stride_pattern = (AccessPattern*)malloc(sizeof(AccessPattern)) ) == NULL){
					printf("Failed to allocate space for Access Pattern!");
					return -1;
				}

				stride_pattern->k = 1;
				stride_pattern->reqSize = list_i->size;
				stride_pattern->startPos = list_i->offset;
				stride_pattern->endPos = list_i->offset;
				stride_pattern->operation = list_i->op;
				stride_pattern->patternType = KD_STRIDED;
				stride_pattern->strideSize[0] = tmp_stride_size;
				stride_pattern->recordNum[0] = 1;
				stride_pattern->mpiRank = mpirank;
				stride_pattern->startTime = list_i->startTime;
				stride_pattern->endTime = list_i->endTime;

				DL_APPEND(*pattern_head, stride_pattern);

				// Delete the records corresponding to this pattern, feed the rest to this pattern
				TraceList *old = list_i;
				TraceList *oldnext = old->next;
				DL_DELETE(*tracelist, old);

				old = oldnext;
				oldnext = old->next;
				while(1){
					if(stride_check(old, stride_pattern)){
						// old in the pattern
						stride_pattern->endPos = old->offset;
						//stride_pattern->endTime = old->op_time;
						stride_pattern->recordNum[0] ++;
						DL_DELETE(*tracelist, old);

					}
					old = oldnext;
					if(old == NULL)
						break;
					oldnext = old->next;
				}

				// add to look up hash table
				UT_lookup *tmp = (UT_lookup*)malloc( sizeof(UT_lookup) );
			    memset(tmp, 0, sizeof(UT_lookup));
				tmp->key.op = stride_pattern->operation;
				tmp->key.mpirank = stride_pattern->mpiRank;
				tmp->key.off = stride_pattern->endPos + stride_pattern->strideSize[0];
				tmp->pattern = stride_pattern;

			    HASH_ADD(hh, lookup, key, sizeof(Lookup_key), tmp);

				// if one stride pattern is found then return
				return stride_pattern->recordNum[0];
			}

			list_j = list_j->next;

		}

		list_i = list_i->next;
	}

	return 0;
}


int check_pattern_same(AccessPattern *p1, AccessPattern *p2)
{
	int i;
	for(i = 0; i < p1->k ; i++){
		if(p1->patternType != p2->patternType || p1->strideSize[i] != p2->strideSize[i] ||p1->mpiRank != p2->mpiRank
				|| p1->k != p2->k)
			return 0;

	}

	return 1;
}

int merge_kd(AccessPattern **kdpattern_head, int k)
{

	AccessPattern *pattern_i = *kdpattern_head;
	AccessPattern *pattern_j = NULL;
	AccessPattern *pattern_k = NULL;
	AccessPattern *pattern_t = NULL;

	int tmpstridesize = 0;
	int tmpstartpos = 0;
	int hasmerge = 0;
	int mergecnt = 0;

	while(pattern_i != NULL){

		//for each stride pattern, try to fit other stride patterns with same stride size to form kd-stride pattern
		if(pattern_i->patternType == KD_STRIDED  && pattern_i->k == k - 1 ){
			pattern_j = pattern_i->next;

			while(pattern_j != NULL){
				hasmerge = 0;
				int checkflag = 0;
				int i;
				for(i = 0; i < pattern_i->k; i++){
					//check pattern_i&j's 1&2-D corresponding value match
					if(!check_pattern_same(pattern_i,pattern_j) || pattern_j->reqSize != pattern_i->reqSize
							|| pattern_i->recordNum[i] != pattern_j->recordNum[i]){
						//
						checkflag = 1;
						break;
					}
				}

				if(checkflag == 1){
					pattern_j = pattern_j->next;
					continue;
				}
				else{
					tmpstridesize = pattern_j->startPos - pattern_i->startPos;
					pattern_k = pattern_j->next;
					while(pattern_k != NULL){
						checkflag = 0;
						int i;
						if(pattern_j == NULL)
							break;
						for(i = 0; i < pattern_j->k; i++){
							//check pattern_j&k's 1&2-D corresponding value match
							if(!check_pattern_same(pattern_j,pattern_k) || pattern_k->reqSize != pattern_j->reqSize
									|| pattern_j->recordNum[i] != pattern_k->recordNum[i]){
								checkflag = 1;
								break;
							}
						}

						if(checkflag == 1){
							pattern_k = pattern_k->next;
							checkflag = 0;
							continue;
						}
						else{
							//found 3 (k-1)D-strided pattern
							if(hasmerge == 0){
								//upgrade pattern_i to 3D pattern and remove pattern_j
								pattern_i->k++;
								pattern_i->recordNum[pattern_i->k-1] = 3;
								pattern_i->strideSize[pattern_i->k-1] = pattern_j->startPos - pattern_i->startPos;
								pattern_i->endPos = pattern_k->endPos;
								tmpstartpos = pattern_k->startPos;
								hasmerge = 1;
								mergecnt++;
							}
							else{
								if(pattern_k->startPos - tmpstartpos != pattern_i->strideSize[pattern_j->k]){
									pattern_j = pattern_j->next;
									continue;
								}
								pattern_i->recordNum[pattern_i->k-1]++;
								pattern_i->endPos = pattern_k->endPos;
								tmpstartpos = pattern_k->startPos;
							}

							hasmerge = 1;
						}
						if(hasmerge!=0)
							pattern_j = pattern_k;
						pattern_k = pattern_k->next;
					}

				}

				//delete patterns that were merged
				AccessPattern *pattern_tmp = NULL;
				int tmpstartpos = pattern_i->startPos;
				int tmp_rank = pattern_i->mpiRank;

				if(hasmerge!=0){
					DL_FOREACH_SAFE(*kdpattern_head, pattern_t, pattern_tmp){
						if(pattern_t->mpiRank == tmp_rank && pattern_t->patternType == KD_STRIDED
										&& pattern_t->startPos == tmpstartpos + tmpstridesize){
							DL_DELETE(*kdpattern_head, pattern_t);
							tmpstartpos = tmpstartpos + tmpstridesize;
						}
					}

				}
				else
					pattern_j = pattern_j->next;

			}

		}
		if(pattern_i == NULL)
			break;
		pattern_i = pattern_i->next;

	}//while pattern_i

	return mergecnt;
}




