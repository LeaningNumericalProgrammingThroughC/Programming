#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define N 17     //標本点数

double f(double x) {
	return 1.0 / (1.0 + 25.0 * x * x);
}

// ラグランジュ補間---------------------------------------------------------
double lagrange(double sample_x[], double sample_y[], int n, double x) {

	double s = 0.0, p;
	int i, j;

	// Σ f(x)*l(x)のループ-------------------------------------------
	for (i = 0; i < n; i++) {

		// Π(x-x_j)/(x_i-x_j)のループ--------------------------------
		p = sample_y[i];
		for (j = 0; j < i; j++) {
			p *= (x - sample_x[j]) / (sample_x[i] - sample_x[j]);
		}
		for (j = i + 1; j < n; j++) {
			p *= (x - sample_x[j]) / (sample_x[i] - sample_x[j]);
		}
		s += p;
	}

	return s;
}

// メイン関数---------------------------------------------------------------
int main(void) {

	double sample_x[N], sample_y[N], x, a = -1.0, b = 1.0;
	int i, m = 100;
	double h = (b - a) / m;

	// 標本点をセット
	for (int i = 0; i < N; i++) {
		sample_x[i] = cos((M_PI * (2. * (i + 1.) - 1.)) / (2. * N));
		sample_y[i] = f(sample_x[i]);
	}

	// 補間と結果出力------------------------------------------------
	printf("    x               f(x)           P_n(x) \n");
	printf("--------------------------------------------\n");
	x = a;
	for (i = 0; i <= m; i++) {

		printf("%10.6lf      %10.6lf      %10.6lf\n", x, f(x), lagrange(sample_x, sample_y, N, x));
		x += h;
	}

	return 0;
}