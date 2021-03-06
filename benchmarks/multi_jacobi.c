#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

#define PRINTMAT 0
#define DEBUG 0

void print_x(int iter, int n, double* x)
{
    int i;
    // append results
    FILE *fp = fopen("x.dat", "a+");
    fprintf(fp, "itereation: %d, ", iter);
    fprintf(fp, "x=[");
    for(i= 0; i < n; i++) {
        fprintf(fp, " %lf", x[i]);
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

    double compute_time, compute_time_start, compute_time_end;
    double io_time, io_time_start, io_time_end;
    double total_time_start, total_time_end;
    double sum, temp, diff, bb, allbb;
    double e = 0.00001;
    double *x, *myx;

    // whether to use prefetch or not ,see below
    int compare_mode = 0;

    // Init MPI
    MPI_Init(&argc, &argv);

    // record start time of program
    total_time_start = MPI_Wtime();

    // get the number of procs and rank in the comm
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if(argc < 4) {
        printf("Usage: jacobi mat.bin dim mat_num [enable prefetch]\n");
        return -1;
    }
    else if(argc == 5) {
        // compare mode
        // 0 means no prefetch, 1 means prefetch
        compare_mode = atoi(argv[4]);
        if(my_rank == 0 ) {
            if(compare_mode == 1)
                printf("Using prefetching\n");
            else
                printf("No prefetching\n");
        }
    }

    // broadcast prefetch mode
    MPI_Bcast( &compare_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // n is dimension of the input matrix
    int n = atoi(argv[2]);

    // each proc get myrows rows of the matrix
    int myrows = n / proc_num;

    // CACHE Allocation
    int *prefetch_cache = malloc( myrows * (n+1) * sizeof(int));

    // allmat is 1-D array to store matrix
    allmat = (int*) malloc( myrows * (n+1) * sizeof(int));

    // mymat makes it a 2-D array
    mymat = (int**) malloc(myrows * sizeof(int *));
    for (i = 0; i < myrows; i++) {
        mymat[i]= &allmat[i * (n+1)];
    }

    // x stores global result, myx stores local result
    x = (double*) malloc( n * sizeof(double) );
    myx = (double*) malloc( myrows * sizeof(double));

    // how many matrices to solve in total
    int mat_num = atoi(argv[3]);

#if DEBUG
    printf("I'm proc%d, n=%d, myrows=%d, mat_num=%d\n", my_rank, n, myrows, mat_num);
    if(my_rank == -1) {
        int dowait = 1;
        while(dowait) {
            ;
        }
    }
#endif

    // File opening
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if(fh==NULL) {
        printf("File not exist\n");
        return -1;
    }

    io_time = 0.0;
    compute_time = 0.0;
    double start, finish;
    // each round read entire matrix and prefetch next round
    for(k = 0; k < mat_num; k++) {

        MPI_Barrier(MPI_COMM_WORLD);
        
        // I/O time
        io_time_start = MPI_Wtime();

        if(compare_mode == 1) {
            // use prefetch

            if(k == 0) {
                start = MPI_Wtime();
                // first time read, no pattern information known so just normal read
                MPI_File_read_at( fh, (myrows * my_rank + k * n) * (n+1) * sizeof(int)
                                  , allmat, myrows * (n+1), MPI_INT, &status );
                finish  = MPI_Wtime();
                if(my_rank ==0) printf("First read time %lf\n",finish - start);

                // According to previous read, predict the next read
                // use non-blocking read so computation could be performed at the same time
                MPI_File_iread_at( fh, (myrows * my_rank + (k+1) * n) * (n+1) * sizeof(int)
                                   , prefetch_cache, myrows * (n+1), MPI_INT, &request );
            }
            else {
                start = MPI_Wtime(); 

                // wait for previous iread(prefetched the predicted next access) completed
                MPI_Wait(&request, &status);
                finish  = MPI_Wtime();
                if(my_rank ==0) printf("Wait time %lf\n",finish - start);

                start = MPI_Wtime();
                // copy prefetched data from cache to target
                memcpy(allmat, prefetch_cache, myrows * (n+1) * sizeof(int));
                finish  = MPI_Wtime();
                if(my_rank ==0) printf("Memcpy time %lf\n",finish - start);


                // next read is predicted so perform prefetch
                if(k != mat_num -1) {
                    MPI_File_iread_at( fh, (myrows * my_rank + (k+1) * n) * (n+1) * sizeof(int)
                                       , prefetch_cache, myrows * (n+1), MPI_INT, &request );
                }
            }

        }
        else {
            // normal read
            MPI_File_read_at( fh, (myrows * my_rank + k * n) * (n+1) * sizeof(int)
                              , allmat, myrows * (n+1), MPI_INT, &status );
        }

        MPI_Barrier(MPI_COMM_WORLD);
        if(my_rank == 0) {
            io_time_end = MPI_Wtime();
            printf("I/O time of %d round: %lf\n",k , io_time_end - io_time_start);
            io_time += io_time_end - io_time_start;
        }

#if PRINTMAT
        // print matrix
        printf("rank %d:\n",my_rank);
        for(i=0; i<myrows; i++) {
            for(j=0; j<n+1; j++) {
                printf(" %4d", mymat[i][j]);
            }
            printf("\n");
        }
#endif


        // set local and global x to zero for each iteration
        memset( myx, 0, sizeof(myx[0]) * myrows );
        memset( x, 0, sizeof(x[0]) * myrows );

        compute_time_start = MPI_Wtime();
        // start  iteration of computation till converge
        iter=0;
        do {
            bb=0.0;
            // all proc get all x
            MPI_Allgather(myx, myrows, MPI_DOUBLE, x, myrows, MPI_DOUBLE, MPI_COMM_WORLD);
            for(i = 0; i < myrows; i++) {

                sum=0.0;
                for(j = 0; j < n; j++) {
                    if(j!=i+myrows*my_rank) {
                        sum=sum+(double)mymat[i][j]*x[j];
                    }
                }

                temp=( (double)mymat[i][n]-sum ) / (double)mymat[i][i+myrows*my_rank];

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

        } while(allbb>=e);

        // gather final x for print of each matrix
        MPI_Allgather(myx, myrows, MPI_DOUBLE, x, myrows, MPI_DOUBLE, MPI_COMM_WORLD);

        if(my_rank ==0 ) {
            // record end time of computation
            compute_time_end = MPI_Wtime();
            printf("Compute time of %d round: %lf\n",k , compute_time_end - compute_time_start);
            compute_time += compute_time_end - compute_time_start;

            // append result to file
            print_x(iter, n, x);
        }


    }//k

    MPI_File_close( &fh );

    if(my_rank == 0) {
        printf("Total I/O time: %lf  bandwidth: %.2lf MB/s\n", io_time, n*(n+1)*4.0/(io_time*1024*1024));
        printf("Total compute time: %lf\n", compute_time);

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

