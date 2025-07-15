#include <stdio.h>

#define N 4

void transpose_matrix(double a[][N]) {

	double tmp;
	int i, j;

	for (i = 0; i < N; i++) {

		for (j = 0; j < i; j++) {
			tmp = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = tmp;
		}
	}
}

int main() {

	double a[N][N] = {
		{ 1.0,  2.0,  3.0,  4.0},
		{ 5.0,  6.0,  7.0,  8.0},
		{ 9.0, 10.0, 11.0, 12.0},
		{13.0, 14.0, 15.0, 16.0}
	};
	int i, j;

	transpose_matrix(a);

	printf("“]’us—ñF\n");
	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			printf("%4.1lf ",a[i][j]);
		}
		printf("\n");
	}

	return 0;
}