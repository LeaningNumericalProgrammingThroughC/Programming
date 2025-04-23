#include <stdio.h>

#define N 11    // 標本点数

// スプライン補間-----------------------------------------------------------
double spline(double sample_x[], double sample_y[], int n, double x) {

	double h[N], w;
	double a[N], b[N], c[N+1], d[N];
	double alpha[N+1], beta[N], gamma[N+1], delta[N+1];

	int i, j, k;

	// 各区間幅h[i]の計算-------------------------------------------------
	for (i = 1; i < N; i++) {
		h[i] = sample_x[i] - sample_x[i - 1];
	}

	// c[i]を求めるための係数行列の設定-----------------------------------
	// 1行目の係数をセット
	alpha[1] = 2 * h[1];
	beta[1] = h[1];
	delta[1] = 3 * (sample_y[1] - sample_y[0]);

	// 2～n行目の係数をセット
	for (i = 2; i < N; i++) {
		alpha[i] = 2 * (h[i - 1] + h[i]);
		beta[i] = h[i - 1];
		gamma[i] = h[i];
		delta[i] = 3 * (h[i] / h[i-1] * (sample_y[i-1] - sample_y[i-2])
			           + h[i-1] / h[i] * (sample_y[i] - sample_y[i-1]));
	}

	// n+1行目の係数をセット
	alpha[N] = 2 * h[N-1];
	gamma[N] = h[N-1];
	delta[N] = 3 * (sample_y[N-1] - sample_y[N-2]);

	// 前進消去----------------------------------------------------------
	for (i = 1; i < N; i++) {
		w = gamma[i + 1] / alpha[i];
		alpha[i + 1] -= w * beta[i];
		delta[i + 1] -= w * delta[i];
	}

	// 後退代入----------------------------------------------------------
	c[N] = delta[N] / alpha[N];
	for (i = N - 1; i >= 1; i--) {
		c[i] = (delta[i] - beta[i] * c[i + 1]) / alpha[i];
	}

	// a,b,dの各係数を計算-----------------------------------------------
	for (i = 1; i < N; i++) {
		a[i] = (1.0 / (h[i]*h[i])) * (c[i] + c[i+1] 
			      - 2.0 * (sample_y[i] - sample_y[i-1]) / h[i]);
		b[i] = (1.0 / h[i]) * (-2.0*c[i] - c[i+1]
			      + 3.0 * (sample_y[i] - sample_y[i-1]) / h[i]);
		d[i] = sample_y[i-1];
	}

	// xがどの区間に含まれるかチェック-----------------------------------
	j = 1;
	while (sample_x[j] < x && j < N-1) j++;

	// S(x)を計算--------------------------------------------------------
	w = x - sample_x[j-1];
	return  a[j] * w * w * w + b[j] * w * w + c[j] * w + d[j];
}

// メイン関数---------------------------------------------------------------
int main() {

	// 標本点の定義
	double sample_x[N] = { 1.0, 2.7, 3.0, 4.5, 5.2,
		                 6.1, 7.4, 8.2, 9.3, 10.0,11.5 };
	double sample_y[N] = { 3.5, 5.9, 6.5, 5.0, 4.1,
		                 3.2, 3.0, 4.3, 9.5, 8.3, 8.2 };
	int n = 100;
	double x, y, h = (sample_x[N - 1] - sample_x[0]) / n;

	// 結果出力----------------------------------------------------------
	printf("    x               y\n");
	printf("----------------------------------\n");
	for (int i = 0; i <= n; i++) {
		x = sample_x[0] + i * h;
		y = spline(sample_x, sample_y, N, x);
		printf("%10.6lf\t%10.6lf\n", x, y);
	}

	return 0;
}