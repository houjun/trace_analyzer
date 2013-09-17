#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#define SIZE 268435456
//#define SIZE 1073741824
int main (int argc, char *argv[])
{
    
    int proc_num, my_rank, len;
    char  hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Status status;
    MPI_Request request;
    MPI_File fh;
    int i, j;

    MPI_Init(&argc, &argv);
    
    // get the number of procs and rank in the comm
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    
    // get host name
    //MPI_Get_processor_name(hostname, &len);
    
    // flush_size in GB
    int flush_size = atoi(argv[1]);

    printf("Will allocate %d GB in total\n", flush_size);
    flush_size = flush_size * 1024.0*1024.0*1024.0 / SIZE;

    char **mem = (char**)malloc(flush_size*sizeof(char*));
    for(i = 0; i < flush_size; i++)
        mem[i] = malloc(SIZE*sizeof(char));

    // each node allocates 1G memory until use up memory.
    for(i = 0; i < flush_size; i++){
        mem = malloc(SIZE*sizeof(char));
        if(mem==NULL){
            printf("allocated %.2lf GB, no more available space\n", (double)(i * (SIZE/(1024.0*1024.0*1024.0))) );
            break;
        }
        memset(mem, 0, SIZE);
    }


    char tmp;
    printf("Press any key to terminate the program\n");
    tmp = getchar();
    

    for(i = 0; i < flush_size; i++)
        free(mem[i]);
    free(mem);

    MPI_Finalize();

    return 0;
}

