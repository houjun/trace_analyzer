/*
 ============================================================================
 Name        : TraceAnalyzer.c
 Author      : Houjun Tang
 Version     :
 Copyright   : Your copyright notice
 Description : Trace Analyzer in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "patterndef.h"
#include "utlist.h"
#include "uthash.h"


int qsort_compare (const void *a, const void *b) {
  return ( *(int*)a - *(int*)b );
}

int main(int argc, char *argv[]) {

	char pname[10][20] = {"CONTIGUOUS","REPEAT","KD_STRIDED","SEQUENTIAL"};
	char oname[4][32] = {"ADIO_READ","ADIO_WRITE","ADIO_READS","ADIO_WRITES"};
	int totalbytes = 0;
	int totalpatterncnt = 0;
	int i = 0;

	// store all traces
	TraceList *adio_write_list_head = NULL;
	TraceList *adio_read_list_head = NULL;
	AccessPattern *pattern_head = NULL, *list_i;

	// look up hash table
	UT_lookup *lookup = NULL;

	if(argc != 3){
		printf("Usage: analyzer input_filename output_filename\n");
		return -1;
	}

	//read from file and perform analysis
	read_radar(argv[1], &adio_write_list_head, &adio_read_list_head, &pattern_head, &totalbytes, &lookup);

	//merge kd patterns
	for(i = 2; i <= PATTERN_K_SIZE_MAX; i++){
		if(merge_kd(&pattern_head,i) == 0)
			break;
	}


	// calculate and output total pattern number
	DL_FOREACH(pattern_head, list_i)
		totalpatterncnt++;

	FILE *fp = NULL;
	if( (fp = fopen(argv[2], "w+")) == NULL) {
		printf("No such file\n");
		return -1;
	}
	if (fp == NULL)
		printf("Error Reading File\n");

	//output pattern list
	fprintf(fp,"%d\n%d\nOP            Pattern \t   MPIRank  StartTime    EndTime    #      Start    End   AccessSize #records StrideSize\n", totalbytes, totalpatterncnt);
	DL_FOREACH(pattern_head, list_i){
		if(list_i->patternType == SEQUENTIAL){
			float five_num_summary[5];
			qsort (list_i->reqOffesets, list_i->recordNum[0], sizeof(int), qsort_compare);
			int req_arr_size = list_i->recordNum[0]-1;
			five_num_summary[0] = list_i->reqOffesets[0];
			five_num_summary[4] = list_i->reqOffesets[req_arr_size];
			if(list_i->recordNum[0] % 2 == 0){//even
				five_num_summary[2] = (list_i->reqOffesets[req_arr_size/2] + list_i->reqOffesets[req_arr_size/2 + 1])/2;
				five_num_summary[1] = list_i->reqOffesets[req_arr_size/4];
				five_num_summary[3] = list_i->reqOffesets[req_arr_size*3/4 + 1];
			}
			else{
				five_num_summary[2] = list_i->reqOffesets[req_arr_size/2];
				five_num_summary[1] = (list_i->reqOffesets[req_arr_size/4] + list_i->reqOffesets[req_arr_size/4 + 1]) / 2;
				five_num_summary[3] = (list_i->reqOffesets[req_arr_size*3/4] + list_i->reqOffesets[req_arr_size*3/4 + 1]) / 2;
			}

			fprintf(fp,"%10s  %12s %4d %12lf %12lf  1 %8d %8d %8d %8d %8d"
										,oname[list_i->operation - 1], pname[list_i->patternType - 1], list_i->mpiRank, list_i->startTime, list_i->endTime
										, list_i->startPos, list_i->endPos,	0, list_i->recordNum[0] ,list_i->strideSize[0]);
			fprintf(fp," \t (%f, %f, %f, %f, %f) \n"
					,five_num_summary[0],five_num_summary[1],five_num_summary[2],five_num_summary[3],five_num_summary[4]);
		}
		// for strided
		else if(list_i->k < 2){
			fprintf(fp,"%10s  %12s %4d %12lf %12lf  1 %8d %8d %8d %8d %8d\n"
							,oname[list_i->operation - 1], pname[list_i->patternType - 1], list_i->mpiRank, list_i->startTime, list_i->endTime, list_i->startPos, list_i->endPos,
							list_i->reqSize, list_i->recordNum[0] ,list_i->strideSize[0]);
		}
		else{
			fprintf(fp,"%10s  %12s\t%4d %12lf %12lf   %d %8d %8d\t",oname[list_i->operation - 1],"KD_STRIDED", list_i->mpiRank, list_i->startTime, list_i->endTime, list_i->k
							, list_i->startPos, list_i->endPos);

			for(i=list_i->k-1; i >= 0;i--){
				fprintf(fp,"(%d, %d, %d)",list_i->reqSize,list_i->recordNum[i],list_i->strideSize[i]);
			}
			fprintf(fp,"\n");
		}
	}



	//clean up to analyze next process's trace records
	AccessPattern *elt, *etmp;
	DL_FOREACH_SAFE(pattern_head, elt, etmp) {
		  DL_DELETE(pattern_head, elt);
	}


	fclose(fp);

	return EXIT_SUCCESS;
}
