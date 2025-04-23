#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#define DIM_Y 2                       // yの次元数
#define DIM_Z (DIM_Y + DIM_Y*DIM_Y)   // ヤコビ行列も含めた次元数
#define N_DIV 100                     // ルンゲ・クッタ法の分割数
#define EPS 1.0e-8                    // ニュートン法の収束判定用閾値
#define MAX_ITER  50                  // ニュートン法の最大反復回数

#define N_DIM 6        // 次元数

// 最大値ノルム------------------------------------------------
double vector_norm_max(double x[]) {

	double max;
	int i;

	max = fabs(x[0]);
	for (i = 1; i < DIM_Y; i++) {
		if (fabs(x[i]) > max) {
			max = fabs(x[i]);
		}
	}

	return max;
}

// 部分ピボット選択付きガウスの消去法----------------------------------------------
void pivoted_gauss(double a[][DIM_Y], double b[], double x[])
{
	int i, j, k, max_i;
	double m, s, max, tmp;

	// 前進消去------------------------------------------------------------
	for (i = 0; i < DIM_Y - 1; i++)
	{
		// 部分ピボット選択-----------------------------------------
		max = fabs(a[i][i]);
		max_i = i;
		for (j = i + 1; j < DIM_Y; j++) {

			if (fabs(a[j][i]) > max) {
				max = fabs(a[j][i]);
				max_i = j;
			}
		}

		if (max < EPS) {

			printf("係数行列が正則ではありません");
			exit(1);
		} else if (max_i != i) {

			for (j = i; j < DIM_Y; j++) {
				tmp = a[i][j]; a[i][j] = a[max_i][j]; a[max_i][j] = tmp;
			}
			tmp = b[i]; b[i] = b[max_i]; b[max_i] = tmp;
		}

		// 前進消去-------------------------------------------------
		for (j = i + 1; j < DIM_Y; j++) {

			m = a[j][i] / a[i][i];

			for (k = i + 1; k < DIM_Y; k++) {
				a[j][k] = a[j][k] - m * a[i][k];
			}
			b[j] = b[j] - m * b[i];
		}
	}

	// 後退代入------------------------------------------------------------
	x[DIM_Y - 1] = b[DIM_Y - 1] / a[DIM_Y - 1][DIM_Y - 1];
	for (i = DIM_Y - 2; i >= 0; i--) {

		s = 0;
		for (j = i + 1; j < DIM_Y; j++) {
			s += a[i][j] * x[j];
		}
		x[i] = (b[i] - s) / a[i][i];
	}
}

// 多変数ルンゲ・クッタ法------------------------------------------------
void runge_kutta_system(double a, double b, int n, double y[],
						       void (*f)(double, double[], double[])) {
	double x = a, h = (b - a) / n;
	double k1[N_DIM], k2[N_DIM], k3[N_DIM], k4[N_DIM], tmp[N_DIM];
	int i, j;

	// 見出しと1行目を出力-----------------------------------------
	/*
	printf("    x         ");
	for (j = 0; j < N_DIM; j++) printf("y[%d]        ", j);
	printf("\n");

	printf("%7.3lf   ", x);
	for (j = 0; j < N_DIM; j++) printf("%10.5lf  ", y[j]);
	printf("\n");
	*/

	for (i = 0; i < n; i++) {

		// k1 を計算 -----------------------------------------
		f(x, y, k1);

		// k2 を計算 -----------------------------------------
		for (j = 0; j < N_DIM; j++) tmp[j] = y[j] + h * k1[j] / 2.0;
		f(x + h / 2.0, tmp, k2);

		// k3 を計算 -----------------------------------------
		for (j = 0; j < N_DIM; j++) tmp[j] = y[j] + h * k2[j] / 2.0;
		f(x + h / 2.0, tmp, k3);

		// k4 を計算 -----------------------------------------
		for (j = 0; j < N_DIM; j++) tmp[j] = y[j] + h * k3[j];
		f(x + h, tmp, k4);

		// y を計算	------------------------------------------
		for (j = 0; j < N_DIM; j++)
			y[j] = y[j] + h / 6.0
			* (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);

		x += h;

		// 結果出力	------------------------------------------
		/*
		printf("%7.3lf   ", x);
		for (j = 0; j < N_DIM; j++) printf("%10.5lf  ", y[j]);
		printf("\n");
		*/
	}
}

// 微分方程式の右辺を定義---------------------------------------------------------
void f(double x, double y[], double r[]) {

	r[0] = y[1];
	r[1] = -0.1 * y[1] - y[0] * y[0] * y[0] + 0.3 * cos(x);
	r[2] = y[3];
	r[3] = -3.0 * y[0] * y[0] * y[2] - 0.1 * y[3];
	r[4] = y[5];
	r[5] = -3.0 * y[0] * y[0] * y[4] - 0.1 * y[5];
}

// シューティング法---------------------------------------------------------------
void shooting(double a, double b, double y_a[], 
	           void (*func)(double, double[], double[]), int* r) {

	double F[DIM_Y], dF[DIM_Y][DIM_Y], dy[DIM_Y], z[DIM_Z];
	int iter = 0;

	// ニュートン法を実行--------------------------------------------------
	while (iter++ < MAX_ITER ) {

		// 各関数の初期値をセット
		z[0] = y_a[0];
		z[1] = y_a[1];
		z[2] = 1; z[3] = 0;
		z[4] = 0; z[5] = 1;

		// ルンゲ・クッタ法でy(b)とヤコビ行列dφ(b,α)/dαを計算
		runge_kutta_system(a, b, N_DIV, z, func);

		// F(α)を計算
		F[0] = z[0] - y_a[0];
		F[1] = z[1] - y_a[1];

		// ヤコビ行列dF/dαを計算
		dF[0][0] = z[2] - 1.0;
		dF[0][1] = z[4];
		dF[1][0] = z[3];
		dF[1][1] = z[5] - 1.0;

		// ニュートン法の公式に従いF(α)の符号を逆にする
		for (int j = 0; j < DIM_Y; j++) F[j] = -F[j];

		// ガウスの消去法でニュートン法の修正量を求める
		pivoted_gauss(dF, F, dy);

		// 近似解の更新
		for (int j = 0; j < DIM_Y; j++) y_a[j] += dy[j];

		// 収束判定
		if (vector_norm_max(dy) <= EPS) break;
	}

	// 戻り値をセット
	*r = iter;
}

// メイン関数---------------------------------------------------------------------
int main() {

	double a = 0.0, b = 2.0 * M_PI, y_a[DIM_Y];
	int iter = 0;

	// 初期値の設定
	y_a[0] = 1.0;
	y_a[1] = 1.0;

	shooting(a, b, y_a, f, &iter);

	// 結果出力------------------------------------------------------------
	if (iter > MAX_ITER) {
		printf("反復回数の上限に達しました．\n");
	} else {
		printf("解の初期値は以下のとおり．（反復回数：%d回）\n", iter);
		for (int j = 0; j < DIM_Y; j++) {
			printf("y_a[%d] = %8.6lf\n", j, y_a[j]);
		}
	}

	return 0;
}