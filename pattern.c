/*
 * pattern.c
 *
 *  Created on: Apr 16, 2013
 *      Author: Houjun
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "patterndef.h"
#include "utlist.h"

int read_radar(char* filename, TraceList** write_list, TraceList** read_list, AccessPattern** access_pattern, int* totalbytesrw){

	char buf[MAX_LINE_LENGTH];
	FILE* fp = NULL;
	//printf("%s\n",filename);
	if( (fp = fopen(filename, "r+")) == NULL) {
	    printf("No such file\n");
	    return -1;
	}
	if (fp == NULL)	{
	    printf("Error Reading File\n");
	}

	int mpirank, filepos, size;
	int i,j, commsize, nrecords, mpairs;
	char operator[MAX_LINE_LENGTH];
	char filepath[MAX_LINE_LENGTH];
	char processtmp[MAX_LINE_LENGTH];
	double timedelta = 0.0;
	int read_count = 0;
	int write_count = 0;
	int opStep = 0;
	int hierarchy_info;


	while(fgets(buf, MAX_LINE_LENGTH, fp)){
		if (*buf == '#') continue; /* ignore comment line */

		if(strcmp(buf,"HEADER\n")==0){
			//read next lines as header
			//<Full File path> <MPI file comm size>
			fgets(buf, MAX_LINE_LENGTH, fp);
			sscanf(buf, "%s %d %d",filepath, &commsize, totalbytesrw);

			while(commsize--){
				//PROCESS <MPI Rank> <n: number of trace entries>
				fgets(buf, MAX_LINE_LENGTH, fp);
				sscanf(buf, "%s %d %d",processtmp, &mpirank, &nrecords);

				for(i=0; i < nrecords; i++){

					//read next n records
					fgets(buf, MAX_LINE_LENGTH, fp);
					//<Op string> <Time delta> <m: # of accesses> [<flag1> <flag2> ... <flagM>]
					//TODO: deal with flags
					sscanf(buf, "%d %s %lf %d",&hierarchy_info, operator, &timedelta, &mpairs);
					for(j=0; j < mpairs; j++){
						fgets(buf, MAX_LINE_LENGTH, fp);
						sscanf(buf, "%d, %d", &filepos, &size);

						//add trace record
						if(strstr(operator,"Read") != NULL){
							TraceList *tmp = addtmp(filepos, size, T_ADIO_READ, opStep++, timedelta, mpirank);
							if(tmp==NULL)
								return -1;
							trace_analyzer(read_list, tmp, access_pattern, &read_count, mpirank);
							//just like below
						}
						else if(strstr(operator,"Write") != NULL){

							TraceList *tmp = addtmp(filepos, size, T_ADIO_WRITE, opStep++, timedelta, mpirank);
							if(tmp==NULL)
								return -1;
							trace_analyzer(write_list, tmp, access_pattern, &write_count, mpirank);
						}

					}

				}

				//perform final analysis when reach end of each process
				trace_analyzer(read_list, NULL, access_pattern, &read_count, mpirank);
				trace_analyzer(write_list, NULL, access_pattern, &write_count, mpirank);


				read_count = 0;
				write_count = 0;

				TraceList* elt;
				TraceList* etmp;
				//clean up to analyze next process's trace records
				DL_FOREACH_SAFE(*read_list,elt,etmp) {
				      DL_DELETE(*read_list,elt);
				}
				DL_FOREACH_SAFE(*write_list,elt,etmp) {
					  DL_DELETE(*write_list,elt);
				}

			}


		}//Header


	}

	fclose(fp);
	return 0;
}

TraceList* addtmp(int filepos, int size, int op, int opStep, double opTime, int mpirank){

	TraceList *tmp ;
	if ( (tmp = (TraceList*)malloc(sizeof(TraceList))) == NULL)
		return NULL;
	tmp->mpirank = mpirank;
	tmp->offset = filepos;
	tmp->size = size;
	tmp->op = op;
	tmp->opStep = opStep;
	tmp->opTime = opTime;
	tmp->next = NULL;
	tmp->prev = NULL;
	return tmp;
}

int repeat_check(TraceList* trace, AccessPattern* pattern){

	if(trace->mpirank == pattern->mpiRank && trace->offset == pattern->startPos && trace->size == pattern->reqSize)
		return 1;
	else
		return 0;
}

int contig_check(TraceList* trace, AccessPattern* pattern){

	//TODO support for same offset but different size
	if(trace->mpirank == pattern->mpiRank && trace->offset == pattern->endPos + pattern->reqSize)
		return 1;
	else
		return 0;
}

int stride_check(TraceList* trace, AccessPattern* pattern){

	//TODO support for same offset but different size
	if(trace->mpirank == pattern->mpiRank && trace->offset == pattern->endPos  + pattern->strideSize[0] && trace->size == pattern->reqSize)
		return 1;
	else
		return 0;
}

int pattern_repeat(TraceList** tracelist, AccessPattern** pattern_head, int mpirank){
	TraceList* list_i = *tracelist;
	TraceList* list_j = NULL;
	TraceList* list_k = NULL;

	int repeatsize = 0;
	int tmpendtime = 0;

	while(list_i != NULL && list_i->next != NULL){
		list_j = list_i->next;
		repeatsize = 0;
		while(list_j != NULL){
			if(list_j->offset == list_i->offset && list_j->size == list_i->size){
				//found one repeat continue to search for more
				tmpendtime = list_j->opStep;
				list_k = list_j->next;
				while(list_k != NULL){
					if(list_k->offset == list_i->offset && list_k->size == list_i->size){
						repeatsize++;
					}

					if(repeatsize >= PATTERN_SIZE_THRESHOLD){
						// add new pattern
						AccessPattern *repeat_pattern;
						if ( (repeat_pattern = (AccessPattern*)malloc(sizeof(AccessPattern)) ) == NULL){
							printf("Failed to allocate space for Access Pattern!");
							return -1;
						}

						repeat_pattern->k = 1;
						repeat_pattern->reqSize = list_i->size;
						repeat_pattern->startPos = list_i->offset;
						repeat_pattern->endPos = list_i->offset;
						repeat_pattern->operation = list_i->op;
						repeat_pattern->patternType = REPEAT;
						repeat_pattern->strideSize[0] = 0;
						repeat_pattern->startStep = list_i->opStep;
						repeat_pattern->recordNum[0] = 1;
						repeat_pattern->endStep = tmpendtime;
						repeat_pattern->mpiRank = mpirank;
						repeat_pattern->startTime = list_i->opTime;

						DL_APPEND(*pattern_head, repeat_pattern);

						// Delete the records corresponding to this pattern, feed the rest to this pattern
						TraceList* old = list_i;
						TraceList* oldnext = old->next;
						DL_DELETE(*tracelist, old);

						old = oldnext;
						oldnext = old->next;
						while(1){
							if(repeat_check(old, repeat_pattern)){
								// old in the pattern
								repeat_pattern->recordNum[0]++;
								//repeat_pattern->endTime = old->op_time;

								DL_DELETE(*tracelist, old);

							}
							old = oldnext;
							if(old == NULL)
								break;
							oldnext = old->next;
						}

						// if one contig pattern is found then return
						return repeat_pattern->recordNum[0];

					}//add pattern

					list_k = list_k->next;
				}

			}

			list_j = list_j->next;
		}

		list_i = list_i->next;
	}

	return 1;
}

int pattern_contig(TraceList** tracelist, AccessPattern** pattern_head, int mpirank){
	TraceList* list_i = *tracelist;
	int tmpendtime = 0;
	int isseq = 0;

	while(list_i != NULL && list_i->next != NULL){

		int contig_size = 1;
		isseq = 0;
		int req_arr_size = 0;

		//check each other records in the list to see if forms contig pattern
		TraceList* list_j = list_i->next;
		int tmppos = list_i->offset;
		int tmpsize = list_i->size;


		while(list_j != NULL){
			if(list_j->offset == tmppos + tmpsize){
				tmpendtime = list_j->opStep;
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
				contig_pattern->startStep = list_i->opStep;
				contig_pattern->endStep = tmpendtime;
				contig_pattern->mpiRank = mpirank;
				contig_pattern->startTime = list_i->opTime;
				contig_pattern->reqOffesets[req_arr_size++] = list_i->size;

				DL_APPEND(*pattern_head, contig_pattern);

				// Delete the records corresponding to this pattern, feed the rest to this pattern
				TraceList* old = list_i;
				TraceList* oldnext = old->next;
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

				// if one contig pattern is found then return
				return contig_pattern->recordNum[0];
			}

		}


		list_i = list_i->next;

	}

	return 0;
}

int pattern_fixed_stride(TraceList** tracelist, AccessPattern** pattern_head, int mpirank){

	TraceList* list_i = *tracelist;
	TraceList* list_j = NULL;
	TraceList* list_k = NULL;
	int tmpendStep = 0;

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
			tmpendStep = list_j->opStep;
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
				stride_pattern->startStep = list_i->opStep;
				stride_pattern->recordNum[0] = 1;
				stride_pattern->endStep = tmpendStep;
				stride_pattern->mpiRank = mpirank;
				stride_pattern->startTime = list_i->opTime;

				DL_APPEND(*pattern_head, stride_pattern);

				// Delete the records corresponding to this pattern, feed the rest to this pattern
				TraceList* old = list_i;
				TraceList* oldnext = old->next;
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

				// if one stride pattern is found then return
				return stride_pattern->recordNum[0];
			}

			list_j = list_j->next;

		}

		list_i = list_i->next;
	}


	return 0;
}

int trace_feed(TraceList** tracelist, TraceList* new_record, AccessPattern** pattern_head, int mpirank){

	// for all known patterns check if new_event could fit in
	// if yes then update the pattern with new_event
	// return (Pattern);
	AccessPattern* pattern = *pattern_head;
	int findpattern = 0;
	while(pattern != NULL){
		if(pattern->patternType == CONTIGUOUS || pattern->patternType == SEQUENTIAL){
			if(contig_check(new_record, pattern)){
				pattern->recordNum[0] ++;
				pattern->endPos = new_record->offset;
				pattern->reqSize = new_record->size;
				//pattern->endTime = new_record->op_time;
				if(pattern->patternType == SEQUENTIAL){
					pattern->reqOffesets[pattern->recordNum[0]-1] = new_record->size;
				}

				findpattern++;
			}
		}
		else if(pattern->patternType == KD_STRIDED){
			if(stride_check(new_record, pattern)){
				pattern->recordNum[0] ++;
				pattern->endPos = new_record->offset;
				//pattern->endTime = new_record->op_time;
				findpattern++;
			}
		}
		else if(pattern->patternType == REPEAT){
			if(repeat_check(new_record, pattern)){
				pattern->recordNum[0] ++;
				pattern->endPos = new_record->offset;
				//pattern->endTime = new_record->op_time;
				findpattern++;
			}
		}
		pattern = pattern->next;
	}

	// if could fit to pattern then don't add to trace list
	if(findpattern == 0){
		if(new_record==NULL){
			printf("Can't feed NULL!\n");
			return 0;
		}

		// if no pattern fit then add to trace list
		DL_APPEND(*tracelist, new_record);

	}


	return findpattern;

}

int merge_markov(AccessPattern** pattern_head, Markovpattern** markov_head){

	AccessPattern* pattern = *pattern_head;
	AccessPattern* pattern_next = NULL;
	Markovpattern* MarkovPattern = NULL;

	while(pattern != NULL){
		int hasmerge = 0;
		MarkovPattern = NULL;

		//for each repeat pattern, try to fit other repeat patterns to form Markov pattern
		if(pattern->patternType == REPEAT){
			pattern_next = pattern->next;

			while(pattern_next != NULL){
				if(pattern_next->patternType != pattern->patternType || pattern_next->strideSize != pattern->strideSize
						|| pattern_next->mpiRank != pattern->mpiRank){
					pattern_next = pattern_next->next;
					continue;
				}
				else{
					if(pattern_next->startStep > pattern->startStep && pattern_next->startStep < pattern->endStep){
						if(MarkovPattern == NULL){
							if ( (MarkovPattern = (Markovpattern*)malloc(sizeof(Markovpattern))) == NULL ){
								printf("Failed to allocate space for Access Pattern!");
								return -1;
							}
							//init
							MarkovPattern->operation = pattern->operation;
							MarkovPattern->startPos[0] = pattern->startPos;
							MarkovPattern->strideSize[0] = pattern->strideSize[0];
							MarkovPattern->reqSize[0] = pattern->reqSize;
							MarkovPattern->recordNum[0] = pattern->recordNum[0];
							MarkovPattern->endPos[0] = pattern->endPos;

							MarkovPattern->startPos[1] = pattern_next->startPos;
							MarkovPattern->strideSize[1] = pattern_next->strideSize[0];
							MarkovPattern->reqSize[1] = pattern_next->reqSize;
							MarkovPattern->recordNum[1] = pattern_next->recordNum[0];
							MarkovPattern->endPos[1] = pattern_next->endPos;

							MarkovPattern->patternType = pattern->patternType;
							MarkovPattern->k = 2;
							MarkovPattern->mpiRank = pattern->mpiRank;
							MarkovPattern->startTime = pattern->startTime;

							DL_APPEND(*markov_head, MarkovPattern);

						}
						else{
							//already allocate that pattern
							MarkovPattern->startPos[MarkovPattern->k] = pattern_next->startPos;
							MarkovPattern->strideSize[MarkovPattern->k] = pattern_next->strideSize[0];
							MarkovPattern->reqSize[MarkovPattern->k] = pattern_next->reqSize;
							MarkovPattern->recordNum[MarkovPattern->k] = pattern_next->recordNum[0];
							MarkovPattern->endPos[MarkovPattern->k] = pattern_next->endPos;

							(MarkovPattern->k)++;
						}

						DL_DELETE(*pattern_head, pattern_next);
						hasmerge = 1;

					}

				}
				pattern_next = pattern_next->next;

			}

			if(hasmerge){
				AccessPattern* tmp = pattern;
				pattern = pattern->next;
				DL_DELETE(*pattern_head, tmp);
			}

		}
		if(pattern == NULL)
			break;
		if(!hasmerge)
			pattern = pattern->next;
	}


	return 1;
}

int check_pattern_same(AccessPattern* p1, AccessPattern* p2){
	int i;
	for(i = 0; i < p1->k ; i++){
		if(p1->patternType != p2->patternType || p1->strideSize[i] != p2->strideSize[i] ||p1->mpiRank != p2->mpiRank || p1->k != p2->k)
			return 0;

	}

	return 1;
}

int merge_kd(AccessPattern** kdpattern_head, int k){

	AccessPattern* pattern_i = *kdpattern_head;
	AccessPattern* pattern_j = NULL;
	AccessPattern* pattern_k = NULL;
	AccessPattern* pattern_t = NULL;

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
				AccessPattern* pattern_tmp = NULL;
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

int trace_analyzer(TraceList** list, TraceList* tmp, AccessPattern** access_pattern, int* count, int mpirank){
	if(tmp != NULL && !trace_feed(list, tmp, access_pattern, mpirank))
		(*count)++;

	if(*count > LIST_SIZE_THRESHOLD || tmp == NULL){
		int decreasecount;
		// perform analysis: call all patterns
		decreasecount = pattern_contig(list, access_pattern, mpirank);
		decreasecount = pattern_fixed_stride(list, access_pattern, mpirank);
		decreasecount = pattern_repeat(list, access_pattern, mpirank);

		*count = *count - decreasecount;
	}

	return 0;
}

