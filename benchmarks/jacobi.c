#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

void print_x(int iter, int n, double* x)
{
    int   i;
    FILE *fp = fopen("x.dat", "w+");
    fprintf(fp, "itereation: %d, ", iter);
    fprintf(fp, "x=[");
    for(i= 0; i < n; i++) {
        fprintf(fp, " %f", x[i]);
    }
    fprintf(fp, " ]");
    fprintf(fp, "\n");
    fclose(fp);
}
int main (int argc, char *argv[])
{

    int proc_num, my_rank;

    MPI_Status status;
    MPI_Request request;
    MPI_File fh;

    int i, j, k, iter;
    int **mymat, *allmat;

    double start_time, end_time, io_time;
    double total_time_start, total_time_end;
    double sum, temp, diff, bb, allbb;
    double e = 0.000001;
    double *x, *myx;

    int split = 4;
    int split_rows;
    int compare_mode = 0;

    // Init MPI
    MPI_Init(&argc, &argv);

    // get the number of procs and rank in the comm
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if(argc < 3){
    	printf("Usage: jacobi mat.bin dim\n");
    	return -1;
    }
    else if(argc == 4){
    	// compare mode
    	// 0 means no prefetch, 1 means prefetch
    	compare_mode = atoi(argv[3]);
    	if(my_rank == 0 ){
    		if(compare_mode == 1)
    			printf("Using prefetching\n");
    		else
    			printf("No prefetching\n");
    	}
    }

    MPI_Bcast( &compare_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // dimension of the input matrix
    int n = atoi(argv[2]);
    int myrows = n / proc_num;
    split_rows = myrows / split;

    // CACHE Allocation
    int *prefetch_cache = malloc( split_rows * (n+1) * sizeof(int));

	if(my_rank == -1) {
		int dowait = 1;
		while(dowait) {
			;
		}
	}


    total_time_start = MPI_Wtime();

    allmat = (int*) malloc( myrows * (n+1) * sizeof(int));
    mymat = (int**) malloc(myrows * sizeof(int *));

    for (i = 0; i<myrows; i++) {
        mymat[i]= &allmat[i * (n+1)];
    }

    x = (double*) malloc((n) * sizeof(double));
    myx = (double*) malloc((myrows) * sizeof(double));
    for(i=0; i<myrows; i++) {
        myx[i]=0.0;
    }

    // split one read to multiple ones
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if(fh==NULL){
    	printf("File not exist\n");
    	return -1;
    }

    io_time = 0.0;
    // each proc read split_rows rows of the matrix, perform first iteration of computation
    for(k = 0; k < split; k++){

        MPI_Barrier(MPI_COMM_WORLD);

        start_time = MPI_Wtime();

        if(compare_mode == 1){
        	// prefetch mode
			if(k == 0){
				// first time read, no pattern information known so just normal read
				MPI_File_read_at(fh, (myrows * my_rank + k * split_rows) * (n+1)  * sizeof(int)
						, &allmat[k * split_rows * (n+1)], split_rows * (n+1), MPI_INT, &status);

				// According to previous read, predict the next read
				// use non-blocking read so computation could be performed at the same time
				MPI_File_iread_at(fh, (n+1)  * sizeof(int) * (myrows * my_rank + (k+1) * split_rows)
									, prefetch_cache, split_rows * (n+1), MPI_INT, &request);
			}
			else{
				// wait for previous iread(prefetched the predicted next access) completed
				MPI_Wait(&request, &status);

				// copy prefetched data from cache to target
				memcpy(&allmat[k * split_rows * (n+1)], prefetch_cache, split_rows * (n+1)* sizeof(int));

				// next read is predicted so perform prefetch
				if(k != split -1){
					MPI_File_iread_at(fh, (n+1)  * sizeof(int) * (myrows * my_rank + (k+1) * split_rows)
									, prefetch_cache, split_rows * (n+1), MPI_INT, &request);
				}
			}
        }
        else{
        	// no prefetch
			MPI_File_read_at(fh, (myrows * my_rank + k * split_rows) * (n+1)  * sizeof(int)
					, &allmat[k * split_rows * (n+1)], split_rows * (n+1), MPI_INT, &status);
        }

        MPI_Barrier(MPI_COMM_WORLD);

		end_time = MPI_Wtime();
		io_time += end_time - start_time;

        bb=0.0;
        // all proc get all x
		MPI_Allgather(myx, myrows, MPI_DOUBLE, x, myrows, MPI_DOUBLE, MPI_COMM_WORLD);
		for(i = k * split_rows; i < k * split_rows + split_rows; i++) {
			sum=0.0;

			for(j=0; j<n; j++) {
				if(j!=i+myrows*my_rank) {
					sum=sum+(double)mymat[i][j]*x[j];
				}
			}

			temp=((double)mymat[i][n]-sum) / (double)mymat[i][i+myrows*my_rank];
			diff=fabs(x[i]-temp);

			if(diff>bb) {
				bb=diff;
			}
			myx[i]=temp;

		}
		// each process get same bb value so all can go out of loop

		MPI_Allreduce( &bb, &allbb, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
    }

    MPI_File_close( &fh );

    if(my_rank == 0)
    	printf("n=%d, reading time: %lf  bandwidth: %.2lf MB/s\n", n, io_time, n*(n+1)*4.0/(io_time*1024*1024));

/*
        printf("rank %d:\n",my_rank);
        for(i=0; i<myrows; i++) {
            for(j=0; j<n+1; j++) {
                printf(" %4d", mymat[i][j]);
            }
            printf("\n");
        }

 */

    // start next iteration of computation till converge
    iter=1;
    double compute_time = MPI_Wtime();
    do {
        bb=0.0;
        // all proc get all x
        MPI_Allgather(myx, myrows, MPI_DOUBLE, x, myrows, MPI_DOUBLE, MPI_COMM_WORLD);
        for(i=0; i<myrows; i++) {
            sum=0.0;

            for(j=0; j<n; j++) {
                if(j!=i+myrows*my_rank) {
                    sum=sum+(double)mymat[i][j]*x[j];
                }
            }

            temp=((double)mymat[i][n]-sum) / (double)mymat[i][i+myrows*my_rank];

            diff=fabs(x[i]-temp);

            if(diff>bb) {
                bb=diff;
            }

            myx[i]=temp;

        }
        // each process get same bb value so all can go out of loop

        MPI_Allreduce( &bb, &allbb, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );

        iter++;
//        if(my_rank == 0)
//            printf("iter = %d, bb = %lf\n",iter, allbb);

    }while(allbb>=e);

    // gather final x for print
    MPI_Allgather(myx, myrows, MPI_DOUBLE, x, myrows, MPI_DOUBLE, MPI_COMM_WORLD);

    if(my_rank ==0 ) {
        // record end time of computation
        end_time = MPI_Wtime();
        printf("Compute time:%lf\n", end_time - compute_time);

        print_x(iter, n, x);
    }

    total_time_end = MPI_Wtime();
    if(my_rank == 0)
    	printf("Total time: %lf\n",total_time_end - total_time_start);

    // free allocated memory

    free(allmat);
    free(mymat);
    free(myx);
    free(x);
    free(prefetch_cache);
    MPI_Finalize();
    return 0;

}

