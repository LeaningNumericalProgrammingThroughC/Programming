#include <stdio.h>

/* 行列の積 */
void matrix_product(double c[][4], double a[][4], double b[][4]) {

	int i, j, k;

	// 行のスープ---------------------------------------------
	for (i = 0; i < 4; i++) {

		// 列のループ-----------------------------------------
		for (j = 0; j < 4; j++) {

			// 内積計算---------------------------------------
			c[i][j] = 0.0;
			for (k = 0; k < 4; k++) {
				c[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}

/* メイン関数 */
int main() {

	double a[4][4] = {
		{ 1.0,  2.0,  3.0,  4.0},
		{ 0.0,  1.0,  2.0, -1.0},
		{-1.0, -3.0,  4.0,  1.0},
		{ 2.0,  0.0,  3.0,  5.0}
	};

	double b[4][4] = {
		{ 4.0, -1.0,  2.0, -2.0},
		{ 1.0,  0.0,  2.0, -3.0},
		{ 2.0,  2.0,  0.0, -1.0},
		{-2.0,  4.0,  1.0, -1.0}
	};

	double c[4][4];
	int i, j;

	matrix_product(c, a, b);

	printf("C=A+Bの計算結果\n");
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			printf("%5.2f,\t", c[i][j]);
		}
		printf("\n");
	}

	return 0;
}