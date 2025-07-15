#include <stdio.h>
#include <math.h>

#define N 4

double frobenius_norm(double a[][N]) {
	
	double s = 0.0;
	int i, j;

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			s += a[i][j] * a[i][j];
		}
	}

	return sqrt(s);
}

int main() {
	double a[N][N] = {
		{ 1.0,  2.0,  3.0,  4.0},
		{ 5.0,  6.0,  7.0,  8.0},
		{ 9.0, 10.0, 11.0, 12.0},
		{13.0, 14.0, 15.0, 16.0}
	};

	printf("フロベニウスノルム：%lf\n", frobenius_norm(a));

	return 0;
}