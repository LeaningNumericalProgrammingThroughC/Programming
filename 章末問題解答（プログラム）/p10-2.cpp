#include <stdio.h>
#include <math.h>
#define N 4

// ベクトルの2ノルム
double vector_norm2(double x[]) {

	double s;
	int i;

	s = 0.0;
	for (i = 0; i < N; i++) {
		s += x[i] * x[i];
	}

	return sqrt(s);
}

// 行列を画面出力
void print_matrix(double a[][N]) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			printf(" %6.3lf ", a[i][j]);

		}
		printf("\n");
	}
}

// QR分解
void qr_decomposition(double q[][N], double r[][N], double a[][N]) {

	double u[N], s;
	int i, j, k;

	// Qの成分を計算する
	for (i = 0; i < N; i++) {

		for (j = 0; j < N; j++) u[j] = a[j][i];

		for (j = 0; j < i; j++) {

			s = 0.0;
			for (k = 0; k < N; k++) {
				s += a[k][i] * q[k][j];
			}
			for (k = 0; k < N; k++) {
				u[k] -= s * q[k][j];
			}
		}

		r[i][i] = vector_norm2(u);
		for (k = 0; k < N; k++) {
			q[k][i] = u[k] / r[i][i];
		}
	}

	// Rの成分を計算する
	for (j = 0; j < N; j++) {

		for (i = 0; i < j; i++) {
			r[i][j] = 0;
			for (k = 0; k < N; k++) {
				r[i][j] += a[k][j] * q[k][i];
			}
		}

		for (i = j + 1; i < N; i++) r[i][j] = 0.0;
	}
}

// メイン関数
int main() {

	double a[N][N] = {
		{ 5.0,  1.0,  3.0,  2.0},
		{ 8.0,  6.0, -2.0,  2.0},
		{ 3.0, -2.0,  7.0,  4.0},
		{-5.0,  3.0, -1.0,  1.0}
	};
	double q[N][N], r[N][N];

	qr_decomposition(q, r, a);

	// 結果出力
	printf("A =\n");
	print_matrix(a);

	printf("\nQ =\n");
	print_matrix(q);

	printf("\nR =\n");
	print_matrix(r);

	return 0;
}