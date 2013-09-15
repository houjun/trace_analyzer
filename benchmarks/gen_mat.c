#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAXVALUE 400

int check_matrix(int** mat,int n)
{
    int i,j;
    int sum;

    for(i=0; i<n; i++) {
        sum=0;
        for(j=0; j<n; j++) {
            if(i!=j)
                sum += fabs(mat[i][j]);
        }

        if(fabs(mat[i][i]) < sum)
            return 0;
    }
    return 1;
}
int main(int argc, char *argv[]) {

        FILE *fout;
        int i,j;

        int n = atoi(argv[1]);

        int *allmat = (int*) malloc( n * (n+1) * sizeof(int));
        int **mat = (int**) malloc(n * sizeof(int *));

        for (i = 0; i<n; i++)
                mat[i]= &allmat[i * (n+1)];

        int iter = 0;
        do{
                iter ++;
                /* Reading the input matrix */
                for(i=0; i<n; i++) {
                        int sum = 0;
                        for(j=0; j<n+1; j++) {
                                mat[i][j] = rand() % MAXVALUE - MAXVALUE/4;
                                if(j != n+1)
                                    sum += fabs(mat[i][j]);
                               // printf("%10d",mat[i][j]);
                        }
                        sum -= mat[i][i];
                        if(sum < 0) sum *= -1;
                        mat[i][i] = sum + 1024 * (rand() % 10);
                        //printf("\n");
                }

                if(iter > 64){
                  printf("Unable to generate disired matrix, exiting...\n");
                  return -1;
                }
        }while(check_matrix(mat,n) != 1);

        printf("\nGot matrix, output to binary file\n");
        fout=fopen("mat.bin", "wb");
        fwrite(allmat, sizeof(allmat[0]), n*(n+1), fout);
        fclose(fout);

        return EXIT_SUCCESS;
}
