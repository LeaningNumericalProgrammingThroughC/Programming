#include <stdio.h>
#include <cmath>

#define N  3            // 次元数
#define EPS 1e-8        // 収束判定用
#define MAX_ITER 100    // 最大反復回数

// 行列とベクトルの積---------------------------------------------------------
void matrix_vector_product(double c[], double a[][N], double b[]) {

	double s;
	int i, j;

	for (i = 0; i < N; i++) {

		s = 0.0;
		for (j = 0; j < N; j++) {

			s += a[i][j] * b[j];
		}
		c[i] = s;
	}
}

// 最大値ノルム----------------------------------------------------------------
double vector_norm_max(double x[]) {
	double max;
	int i;

	max = 0.0;
	for (i = 0; i < N; i++) {
		if (fabs(x[i]) > max) {
			max = fabs(x[i]);
		}
	}

	return max;
}

// LU分解----------------------------------------------------------------------
void lu_decomposition(double a[][N], int p[]) {

	double m, s, max, tmp;
	int i, j, k, max_i;

	for (i = 0; i < N - 1; i++)
	{
		// 部分ピボット選択-----------------------------------
		// 最大行の検索
		max = fabs(a[i][i]);
		max_i = i;
		for (j = i + 1; j < N; j++) {

			if (fabs(a[j][i]) > max) {
				max = fabs(a[j][i]);
				max_i = j;
			}
		}

		// 入れ換えた行を配列に記憶
		p[i] = max_i;

		// 正則性の判定
		if (max < EPS) {
			printf("係数行列が正則ではありません");
			exit(1);
		}

		// 行の入れ換え
		if (max_i != i) {

			for (j = 0; j < N; j++) {
				tmp = a[i][j]; a[i][j] = a[max_i][j]; a[max_i][j] = tmp;
			}
		}

		// 前進消去-------------------------------------------
		for (j = i + 1; j < N; j++) {

			m = a[j][i] / a[i][i];

			// 行列Uの要素を計算
			for (k = i + 1; k < N; k++) {
				a[j][k] = a[j][k] - m * a[i][k];
			}

			// 行列Lの要素を格納
			a[j][i] = m;
		}
	}
}

// 前進代入と後退代入----------------------------------------------------------
void lu_substitution(double a[][N], double b[], int p[]) {

	double tmp, s;
	int i, j, k;

	// 方程式の右辺Pbを求める-----------------------------------
	for (i = 0; i < N - 1; i++) {
		tmp = b[i]; b[i] = b[p[i]]; b[p[i]] = tmp;
	}

	// 前進代入-------------------------------------------------
	for (i = 0; i < N; i++) {

		for (j = 0; j < i; j++) {
			b[i] -= a[i][j] * b[j];
		}
	}

	// 後退代入-------------------------------------------------
	b[N - 1] = b[N - 1] / a[N - 1][N - 1];
	for (i = N - 2; i >= 0; i--) {

		s = 0;
		for (j = i + 1; j < N; j++) {
			s += a[i][j] * b[j];
		}
		b[i] = (b[i] - s) / a[i][i];
	}
}

// LU分解による逆行列の計算-------------------------------------
void lu_inverse_matrix(double c[][N], double a[][N]) {

	double b[N];
	int i, j, p[N];

	// LU分解
	lu_decomposition(a, p);

	for (i = 0; i < N; i++) {

		// ベクトルbを第i成分が1の単位ベクトルに初期化
		for (j = 0; j < N; j++) {
			b[j] = (i == j) ? 1.0 : 0.0;
		}

		// Ax=bを解く．
		lu_substitution(a, b, p);

		// 逆行行列の第i列に解xを代入する．
		for (j = 0; j < N; j++) {
			c[j][i] = b[j];
		}
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

// BU_N型の配列を初期化--------------------------------------------------------
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

// BU_N型の配列から関数値とヤコビ行列を抽出------------------------------------
void split(const BU_N BU_x[], double x[], double d[][N]) {

	for (int i = 0; i < N; i++) {

		x[i] = BU_x[i].v;
		for (int j = 0; j < N; j++) {
			d[i][j] = BU_x[i].d[j];
		}
	}
}

// 関数の定義(BU型用）---------------------------------------------------------
void func_BU(BU_N f[], BU_N x[]) {
	f[0] = 5.0 * x[0] * x[0] * x[0] - 2.0 * x[0] * x[0] + x[1] * x[2] - x[1] - 3.0;
	f[1] = 2.0 * x[1] * x[1] * x[2] + 7.0 * x[0] - 2.0 * x[2] - 8.0;
	f[2] = x[0] * x[1] * x[2] + 5.0 * x[0] * x[1] - 4.0 * x[2] + 1.0;
}

// 関数の定義(double用)--------------------------------------------------------
void func_double(double f[], double x[]) {
	f[0] = 5.0 * x[0] * x[0] * x[0] - 2.0 * x[0] * x[0] + x[1] * x[2] - x[1] - 3.0;
	f[1] = 2.0 * x[1] * x[1] * x[2] + 7.0 * x[0] - 2.0 * x[2] - 8.0;
	f[2] = x[0] * x[1] * x[2] + 5.0 * x[0] * x[1] - 4.0 * x[2] + 1.0;
}

// 簡易ニュートン法（多変数版）------------------------------------------------
void simplified_newton_n(double x[], int* r) {

	BU_N BU_x[N], BU_y[N];
	double f[N], J[N][N], delta_x[N], i_J[N][N];

	// ベクトルxの値を自動微分型変数へ代入
	init(BU_x, x);

	// 自動微分型で関数値と偏微分値を計算
	func_BU(BU_y, BU_x);

	// 関数値と偏微分値を分離して抽出
	split(BU_y, f, J);

	// Jの逆行列を求める
	lu_inverse_matrix(i_J, J);

	int iter = 0;
	do {

		matrix_vector_product(delta_x, i_J, f);

		// 近似解xの更新
		for (int i = 0; i < N; i++) x[i] -= delta_x[i];

		func_double(f, x);

		iter++;
	} while (vector_norm_max(delta_x) > EPS && iter < MAX_ITER);

	// 戻り値をセット
	*r = iter;
}

// メイン関数------------------------------------------------------------------
int main() {

	double x[N] = { 1.0,1.0,2.0 };
	int iter;

	// 簡易ニュートン法
	simplified_newton_n(x, &iter);

	// 結果出力--------------------------------------------------------
	if (iter == MAX_ITER) {
		printf("解が見つかりませんでした．\n");
		return 1;
	}
	else {
		printf("反復回数は%d回，解は以下です．\n", iter);
		for (int i = 0; i < N; i++) {
			printf("x[%d] = %lf\n", i, x[i]);
		}
	}

	return 0;
}