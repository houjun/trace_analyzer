#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#define BLOCK_SZ (1<<26) // 8MB
#define NUM_INTS (BLOCK_SZ/sizeof(int))

static size_t mins(size_t a, size_t b) { return a < b ? a : b; }

void usage(){
    printf( "Usage:\n"
			"gen_4D FILENAME.0~9 B X Y Z t\n");
    exit(1);
}
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 1){
        fprintf(stderr, "run with only one process...\n"); 
        exit(1);
    }
    if (argc < 3){ usage(); }
    
    /* grab required arguments */
    char * fname = argv[1];
    size_t x, y, z, t, b;

    b = atoi(argv[2]);
    x = atoi(argv[3]);
    y = atoi(argv[4]);
    z = atoi(argv[5]);
    t = atoi(argv[6]);

    MPI_Info info = MPI_INFO_NULL;
    MPI_Status status;

    if (argc != 7){
        fprintf(stderr, "%d: Not enough arguments\n", argc);
        usage();
    }
    

    MPI_File fh;
    int err = MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE, info, &fh);
    assert(err == MPI_SUCCESS);

    unsigned int *buf = (unsigned int*)malloc(BLOCK_SZ);
    assert(buf != NULL);

    size_t total = b*x*y*z*t*sizeof(int);
    unsigned int ctr = 0;
    size_t written, to_write;
    for (written = 0; written < total; written+=to_write){
        ctr = 0;
        to_write = mins(BLOCK_SZ, total-written);
        size_t i;
        /* init data in case we want to inspect later */
        for (i = written/sizeof(int); i < to_write/sizeof(int); i++){
            buf[i] = ctr++;
        }
        err = MPI_File_write(fh, buf, (int)to_write, MPI_BYTE, MPI_STATUS_IGNORE);
        assert(err == MPI_SUCCESS);
    }
    err = MPI_File_close(&fh);
    assert(err == MPI_SUCCESS); 

    /*
    // check what has been written
    err = MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    assert(err == MPI_SUCCESS);
    MPI_File_read_at(fh, 0, buf, total, MPI_INT, &status);
    int i;
    for(i = 0; i < b*x*y*z*t; i++){
    	if(i % (x*b) == 0)
    		printf("\n");
    	if(i % (x*b*y) == 0)
    		printf("\n=================\n\n");
    	printf(" %3d",buf[i]);
    }
	printf("\n");

    err = MPI_File_close(&fh);
    assert(err == MPI_SUCCESS);
     */

    MPI_Finalize();
    return 0;
}

