/*
 * read_4D.c
 *
 *  Created on: Sep 25, 2013
 *      Author: houjun
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

void usage(){
    printf( "Usage:\n"
			"read_4D FILENAME B X Y Z T_start T_end T_replay_start T_replay_end\n");
    exit(1);
}

int main(int argc, char *argv[])
{
    int proc_num, my_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    // check arguments
    if (argc != 11){ usage(); }

    int b, x, y, z, t_start, t_end, t_replay_start, t_replay_end, replay_time;
    int i, j, k, t;
    MPI_Status status;

    // init
    char *fname = argv[1];
    b = atoi(argv[2]);				// number of variables(int)
    x = atoi(argv[3]);				// number of rows of cubic
    y = atoi(argv[4]);
    z = atoi(argv[5]);
    t_start = atoi(argv[6]);		// start time step
    t_end = atoi(argv[7]);
    t_replay_start = atoi(argv[8]); // "replay" start time step
    t_replay_end = atoi(argv[9]);
    replay_time = atoi(argv[10]);

    printf("b:%d x:%d y:%d z:%d t_start:%d t_end:%d t_replay_start:%d t_replay_start:%d \n",
    		b, x, y, z, t_start, t_end, t_replay_start, t_replay_end);

    MPI_Info info = MPI_INFO_NULL;
    MPI_File fh;

/*              | b |
 			    ____________
			   /           /|  ...Proc0
			  /           //|
			 /___________// |z
	myrows I |__________|/  |
			 |			|  /
			 |          | /y
			 |          |/    ...ProcN
			 ------------
				  x*b
				  T0
*/

    // distribute work to different procs
    int myrows = y / proc_num;
    int myreadsize = b * x * myrows * z * (t_end - t_start);
    // size of on row
    int myonereadsize = x * b;

    // allocate buffer for reading myrows * x * z * b bytes data
    // for each time step
    int *buf = (int*)malloc(myreadsize * sizeof(int));
    assert(buf != NULL);

    // Open file
    int err = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    assert(err == MPI_SUCCESS);

    int start_offset = 0;
    int read_cnt = 0;
    // First read all from t_start to t_end
    for(t = t_start; t < t_end; t++){
    	// each time step start offset is t*b*x*y*z*sizeof(int)
    	start_offset = t * b * x * y * z;
    	// read a slice of each time step
    	for(i = 0; i < z; i++){
        	for(j = 0; j < myrows; j++){
        		MPI_File_read_at(fh, (start_offset + i*b*x*y + j*b*x + my_rank * myonereadsize * myrows) * sizeof(int)
        				, &buf[read_cnt*myonereadsize], myonereadsize, MPI_INT, &status);
        		read_cnt++;
        	}
    	}

    }

    for(t = t_replay_start; t < t_replay_end; k++){
    	// each time step start offset is t*b*x*y*z*sizeof(int)
		start_offset = t * b * x * y * z;
		// read a slice of each time step
		for(i = 0; i < z; i++){
			for(j = 0; j < myrows; j++){
				MPI_File_read_at(fh, (start_offset + i*b*x*y + j*b*x + my_rank * myonereadsize * myrows) * sizeof(int)
						, &buf[read_cnt*myonereadsize], myonereadsize, MPI_INT, &status);
				read_cnt++;
			}
		}
    }

    err = MPI_File_close(&fh);
    //assert(err == MPI_SUCCESS);

    // check read numbers
    if(my_rank == 1){
		int cnt = 0;
		for(i = 0; i < b*x*myrows*z*(t_end-t_start); i++){
			if(i % (x*b) == 0)
				printf("\n");
			if(i % (x*b*myrows) == 0)
				printf("\n==============%d============\n\n",cnt++);
			printf(" %3d",buf[i]);
		}
		printf("\n");
    }
	free(buf);

	MPI_Finalize();
    return 0;
}


