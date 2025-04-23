#include <stdio.h>
#include <math.h>

#define N_DIM 2        // 次元数
#define N_DIV 100      // 分割数

// 微分方程式の右辺の関数定義--------------------------------------------
void f(double x, double y[], double r[]) {

	r[0] = y[1];
	r[1] = -y[1] + 2 * y[0] + 2 * exp(-x);
}

// 多変数ルンゲ・クッタ法------------------------------------------------
void runge_kutta_system(double a, double b, int n, double y[], 
                         void (*f)(double, double [], double [])) {
	double x = a, h = (b - a) / n;
	double k1[N_DIM], k2[N_DIM], k3[N_DIM], k4[N_DIM], tmp[N_DIM];
	int i, j;

	// 見出しと1行目を出力-----------------------------------------
	printf("    x         ");
	for (j = 0; j < N_DIM; j++) printf("y[%d]        ", j);
	printf("\n");

	printf("%7.3lf   ", x);
	for (j = 0; j < N_DIM; j++) printf("%10.5lf  ", y[j]);
	printf("\n");

	for (i = 0; i < n; i++) {

		// k1 を計算 -----------------------------------------
		f(x, y, k1);

		// k2 を計算 -----------------------------------------
		for (j = 0; j < N_DIM; j++) tmp[j] = y[j] + h*k1[j] / 2.0;
		f(x + h / 2.0, tmp, k2);

		// k3 を計算 -----------------------------------------
		for (j = 0; j < N_DIM; j++) tmp[j] = y[j] + h*k2[j] / 2.0;
		f(x + h / 2.0, tmp, k3);

		// k4 を計算 -----------------------------------------
		for (j = 0; j < N_DIM; j++) tmp[j] = y[j] + h*k3[j];
		f(x + h, tmp, k4);

		// y を計算	------------------------------------------
		for (j = 0; j < N_DIM; j++)
			y[j] = y[j] + h / 6.0
			* (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);

		x += h;

		// 結果出力	------------------------------------------
		printf("%7.3lf   ", x);
		for (j = 0; j < N_DIM; j++) printf("%10.5lf  ", y[j]);
		printf("\n");
	}
}

// メイン関数------------------------------------------------------------
int main() {

	double y_0[N_DIM] = { 1.0, -1.0 };

	runge_kutta_system(0.0, 1.0, N_DIV, y_0, f);

	return 0;
}