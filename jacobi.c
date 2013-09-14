#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

    MPI_Init(&argc, &argv);

    // get the number of procs and rank in the comm
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int i,j,iter;
    int n = 10000;
    int **mymat, *allmat;
    double start_time, end_time, total_time;
    double sum,temp,diff,bb,e;
    int myrows = n / proc_num;

    allmat = (int*) malloc( myrows * (n+1) * sizeof(int));
    mymat = (int**) malloc(myrows * sizeof(int *));

    for (i = 0; i<myrows; i++) {
        mymat[i]= &allmat[i * (n+1)];
    }

    if(my_rank == -1) {
        int dowait = 1;
        while(dowait) {
            ;
        }
    }

    start_time = MPI_Wtime();

//    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, argv[1], MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_set_view( fh, (n+1) * myrows * my_rank * sizeof(int), MPI_INT, MPI_INT, "native", MPI_INFO_NULL );
    MPI_File_read_all(fh, allmat, myrows * (n+1) , MPI_INT, &status);

    MPI_File_close( &fh );

    if(my_rank == 0){
        end_time = MPI_Wtime();
        total_time = end_time - start_time;

        printf("n=%d, reading time:%lf\n", n, total_time);
    }

/*
        printf("rank %d:\n",my_rank);
        for(i=0; i<myrows; i++) {
            for(j=0; j<n+1; j++) {
                printf(" %4d", mymat[i][j]);
            }
            printf("\n");
        }
*/
    
    
    // start computation
    double *x = (double*) malloc((n) * sizeof(double));
    double *myx = (double*) malloc((myrows) * sizeof(double));
    for(i=0; i<myrows; i++) {
        myx[i]=0.0;
    }

    e = 0.0000001;
    iter=0;

    double allbb;

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

    // free allocated memory

    free(allmat);
    free(mymat);
    free(myx);
    free(x);
    MPI_Finalize();
    
    return 0;

}

