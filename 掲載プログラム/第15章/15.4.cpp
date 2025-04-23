#include <stdio.h>
#define _USE_MATH_DEFINES 
#include <cmath>

#define N 2            // 次元数
#define DIV 100        // 分割数
#define EPS 1.0e-8     // 収束判定用閾値
#define MAX_ITER  50   // 最大反復回数

// 最大値ノルム------------------------------------------------
double vector_norm_max(double x[]) {

	double max;
	int i;

	max = fabs(x[0]);
	for (i = 1; i < N; i++) {
		if (fabs(x[i]) > max) {
			max = fabs(x[i]);
		}
	}

	return max;
}

// 部分ピボット選択付きガウスの消去法----------------------------------------------
void pivoted_gauss(double a[][N], double b[], double x[])
{
	int i, j, k, max_i;
	double m, s, max, tmp;

	// 前進消去------------------------------------------------------------
	for (i = 0; i < N - 1; i++)
	{
		// 部分ピボット選択-----------------------------------------
		max = fabs(a[i][i]);
		max_i = i;
		for (j = i + 1; j < N; j++) {

			if (fabs(a[j][i]) > max) {
				max = fabs(a[j][i]);
				max_i = j;
			}
		}

		if (max < EPS) {

			printf("係数行列が正則ではありません");
			exit(1);
		} else if (max_i != i) {

			for (j = i; j < N; j++) {
				tmp = a[i][j]; a[i][j] = a[max_i][j]; a[max_i][j] = tmp;
			}
			tmp = b[i]; b[i] = b[max_i]; b[max_i] = tmp;
		}

		// 前進消去-------------------------------------------------
		for (j = i + 1; j < N; j++) {

			m = a[j][i] / a[i][i];

			for (k = i + 1; k < N; k++) {
				a[j][k] = a[j][k] - m * a[i][k];
			}
			b[j] = b[j] - m * b[i];
		}
	}

	// 後退代入------------------------------------------------------------
	x[N - 1] = b[N - 1] / a[N - 1][N - 1];
	for (i = N - 2; i >= 0; i--) {

		s = 0;
		for (j = i + 1; j < N; j++) {
			s += a[i][j] * x[j];
		}
		x[i] = (b[i] - s) / a[i][i];
	}
}

// 自動微分クラス（多変数版）------------------------------------------------------
class BU_N {
public:

	// データメンバ---------------------------------------------------
	double v;
	double d[N] = { 0.0 };

	// コンストラクタ-------------------------------------------------
	BU_N(double x = 0.0) : v(x) { }

	// 出力関数-------------------------------------------------------
	void print() {
		printf("%lf < ", v);
		for (int i = 0; i < N - 1; i++) {
			printf("%lf, ", d[i]);
		}
		printf("%lf >\n", d[N - 1]);
	}
};

// 自動微分クラスに対する演算子多重定義--------------------------------------------
/* 単項マイナス演算子 */
BU_N operator-(const BU_N& a) {
	BU_N r(-a.v);
	for (int i = 0; i < N; i++) { r.d[i] = -a.d[i]; }
	return r;
}
/* 加算 */
BU_N operator+(const BU_N& a, const BU_N& b) {
	BU_N r(a.v + b.v);
	for (int i = 0; i < N; i++) { r.d[i] = a.d[i] + b.d[i]; }
	return r;
}
/* 減算 */
BU_N operator-(BU_N a, BU_N b) {
	BU_N r(a.v - b.v);
	for (int i = 0; i < N; i++) { r.d[i] = a.d[i] - b.d[i]; }
	return r;
}
/* 乗算 */
BU_N operator*(const BU_N& a, const BU_N& b) {
	BU_N r(a.v * b.v);
	for (int i = 0; i < N; i++) { r.d[i] = a.d[i] * b.v + a.v * b.d[i]; }
	return r;
}
/* 除算 */
BU_N operator/(const BU_N& a, const BU_N& b) {
	BU_N r(a.v / b.v);
	for (int i = 0; i < N; i++) {
		r.d[i] = (a.d[i] * b.v - a.v * b.d[i]) / (b.v * b.v);
	}
	return r;
}

// 初等関数--------------------------------------------------------------------
BU_N sin(const BU_N& x) {
	BU_N r(std::sin(x.v));
	for (int i = 0; i < N; i++) { r.d[i] = std::cos(x.v) * x.d[i]; }
	return r;
}

// BU_N型の配列を初期化------------------------------------------------------------
void init(BU_N BU_x[], const double x[]) {

	for (int i = 0; i < N; i++) {

		/* ベクトルxの値を自動微分クラスの値に代入 */
		BU_x[i].v = x[i];

		/* ヤコビ行列の値を単位行列で初期化 */
		for (int j = 0; j < N; j++) {
			BU_x[i].d[j] = (i == j) ? 1.0 : 0.0;
		}
	}
}

// BU_N型の配列から関数値とヤコビ行列を抽出----------------------------------------
void split(const BU_N BU_x[], double x[], double d[][N]) {

	for (int i = 0; i < N; i++) {

		x[i] = BU_x[i].v;
		for (int j = 0; j < N; j++) {
			d[i][j] = BU_x[i].d[j];
		}
	}
}

// 微分方程式の右辺-------------------------------------------------------------
void func(BU_N f[], double x, BU_N y[]) {

	f[0] = y[1];
	f[1] = -0.1 * y[1] - y[0] * y[0] * y[0] + 0.3 * cos(x);
}

// 自動微分型用のルンゲ・クッタ法-----------------------------------------------
void runge_kutta_system_BU(BU_N f_end[N], BU_N y_0[], double start, double end, 
	void (*func)(BU_N[], double, BU_N[])) {

	BU_N y[N], y_old[N], f[N], y_tmp[N];
	BU_N k1[N], k2[N], k3[N], k4[N];
	double h = (end - start) / DIV, x = start;

	int k = 0;

	for (int i = 0; i < N; i++) {
		y[i] = y_0[i];
	}

	for (k = 0; k < DIV; k++) {

		for (int j = 0; j < N; j++)	y_old[j] = y[j];

		func(f, x, y);
		for (int j = 0; j < N; j++)  k1[j] = f[j];

		for (int j = 0; j < N; j++) y_tmp[j] = y[j] + h / 2. * k1[j];
		func(f, x + h / 2, y_tmp);
		for (int j = 0; j < N; j++) k2[j] = f[j];

		for (int j = 0; j < N; j++) y_tmp[j] = y[j] + h / 2. * k2[j];
		func(f, x + h / 2, y_tmp);
		for (int j = 0; j < N; j++) k3[j] = f[j];

		for (int j = 0; j < N; j++) y_tmp[j] = y[j] + h * k3[j];
		func(f, x + h, y_tmp);
		for (int j = 0; j < N; j++) k4[j] = f[j];

		for (int j = 0; j < N; j++) {
			y[j] = y[j] + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
		}

		x += h;
	}

	for (int j = 0; j < N; j++) f_end[j] = y[j];
}

// シューティング法---------------------------------------------------------------
int shooting_BU(double alpha[], double start, double end,
	void (*func)(BU_N[], double, BU_N[]), int* r) {

	BU_N ad_g[N], ad_alpha[N];
	double g[N], dg[N][N], d_alpha[N];
	int converged = 0, n_iter = 0, i, j;


	while (n_iter < MAX_ITER) {

		// alpha を自動微分型に変換----------------------------------------------
		init(ad_alpha, alpha);

		// ルンゲ・クッタ法でG(α)を計算---------------------------------------------
		runge_kutta_system_BU(ad_g, ad_alpha, start, end, func);
		for (i = 0; i < N; i++) ad_g[i] = ad_g[i] - ad_alpha[i];

		// 自動微分型ad_gを関数値g(α)とヤコビ行列dgに分解-------------------------
		split(ad_g, g, dg);

		// ニュートン法の公式に従いg(α)の符号を逆にする------------------------
		for (j = 0; j < N; j++) g[j] = -g[j];

		// ガウスの消去法により，ニュートン法の修正量を求める--------------------
		pivoted_gauss(dg, g, d_alpha);

		// 近似解の更新---------------------------------------------------------
		for (int j = 0; j < N; j++) alpha[j] += d_alpha[j];

		n_iter++;

		if (vector_norm_max(d_alpha) <= EPS) {
			converged = 1;
			break;
		}
	}

	// 戻り値をセット
	*r = n_iter;
	return converged;
}

// メイン関数---------------------------------------------------------------------
int main() {

	double a = 0.0, b = 2.0 * M_PI, y_a[N];
	int converged, iter = 0;

	// 初期値の設定
	y_a[0] = 1.0;
	y_a[1] = 1.0;

	converged = shooting_BU(y_a, a, b, func, &iter);
	
	// 結果出力------------------------------------------------------------
	if (converged) {
		printf("解の初期値は以下のとおり．（反復回数：%d回）\n", iter);
		for (int j = 0; j < N; j++) {
			printf("y_a[%d] = %8.6lf\n", j, y_a[j]);
		}
	} else {
		printf("反復回数の上限に達しました．\n");
	}

	return 0;
}
