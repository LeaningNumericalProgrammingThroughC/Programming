#include <stdio.h>

#define N_DIM 2        // 次元数
#define N_DIV 100      // 分割数

// 微分方程式の右辺の関数定義--------------------------------------------
void f(double x, double y[], double r[]) {

	r[0] = y[1];
	r[1] = 2.0 * y[0] * y[0] + 2.0 * y[0] - x;
}

// リープフロッグ法------------------------------------------------
void leapfrog_system(double a, double b, int n, double y_i_m_1[],
	void (*f)(double, double[], double[])) {

	double x = a, h = (b - a) / n;
	double f_x[N_DIM], y_i_p_1[N_DIM], y_i[N_DIM];
	int i, j;

	// 見出しと1行目を出力--------------------------------
	printf("    x         ");
	for (j = 0; j < N_DIM; j++) printf("y[%d]        ", j);
	printf("\n-----------------------------------\n");

	printf("%7.3lf   ", x);
	for (j = 0; j < N_DIM; j++) printf("%10.5lf  ", y_i_m_1[j]);
	printf("\n");


	// f_x を計算 -----------------------------------------
	f(x, y_i_m_1, f_x);

	// y を計算	------------------------------------------
	for (i = 0; i < N_DIM; i++) {
		y_i[i] = y_i_m_1[i] + h * f_x[i];
	}

	x += h;

	// 結果出力	------------------------------------------
	printf("%7.3lf   ", x);
	for (j = 0; j < N_DIM; j++) printf("%10.5lf  ", y_i[j]);
	printf("\n");

	for (i = 1; i < n; i++) {

		// y_1,y_2の傾きを計算 -------------------------------
		f(x, y_i, f_x);


		// y を計算	------------------------------------------
		for (j = 0; j < N_DIM; j++) {
			y_i_p_1[j] = y_i_m_1[j] + 2.0 * h * f_x[j];
		}

		x += h;

		// 結果出力	------------------------------------------
		printf("%7.3lf   ", x);
		for (j = 0; j < N_DIM; j++) printf("%10.5lf  ", y_i_p_1[j]);
		printf("\n");

		for (j = 0; j < N_DIM; j++) {
			y_i_m_1[j] = y_i[j];
			y_i[j] = y_i_p_1[j];
		}
	}
}

// メイン関数------------------------------------------------------------
int main() {

	double y_0[N_DIM] = { 1.0, 0.0 };

	// リープ・フロッグ法
	leapfrog_system(0.0, 1.0, N_DIV, y_0, f);

	return 0;
}