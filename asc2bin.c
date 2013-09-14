/*
 ============================================================================
 Name        : test.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */
#include <stdlib.h>
#include <stdio.h>


int main(int argc, char *argv[]) {
	FILE *fin, *fout;
	int i,j,n;
	int **mat, *allmat;

	// read mat file
	fin=fopen(argv[1], "r");

	/* Reading the maxtrix dimension */
	fscanf(fin, "%d",&n);

	// store in row major
	allmat = (int*) malloc( n * (n+1) * sizeof(int));
	mat = (int**) malloc(n * sizeof(int *));

	for (i = 0; i<n; i++)
		mat[i]= &allmat[i * (n+1)];

	/* Reading the input matrix */
	for(i=0; i<n; i++) {
		for(j=0; j<n+1; j++) {
			fscanf(fin, "%d", &mat[i][j]);
		}
	}
	fclose(fin);

	/*
	for(i=0; i<n; i++) {
		for(j=0; j<n+1; j++) {
			printf(" %3lf", mat[i][j]);
		}
		printf("\n");
	}
*/


	fout=fopen("mat.bin", "wb");
	fwrite(allmat, sizeof(allmat[0]), n*(n+1), fout);
	fclose(fout);

	return 0;
}
