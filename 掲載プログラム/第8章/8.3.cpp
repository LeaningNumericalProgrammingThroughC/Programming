#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#define N    3            // 次元数（x:2次元 + t:1次元）
#define DELTA_S    0.1    // 解曲線追跡の刻み幅
#define EPS   1e-10       // ニュートン法の収束判定用
#define MAX_ITER 100      // ニュートン法の最大反復回数
#define MAX_STEP 1000     // ホモトピー法の最大ステップ数

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

// 部分ピボット選択付きガウスの消去法----------------------------------------------
void piboted_gauss(double a[][N], double b[], double x[])
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
		}
		else if (max_i != i) {

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

/* 最大値ノルム */
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

// グローバル変数--------------------------------------------------------------------

// g_f_gamma：f(γ)
double g_f_gamma_1, g_f_gamma_2;
// g_n：解曲線の接ベクトル，かつ平面q(x,t)=0の法線ベクトル
double g_n[N];
// g_p_c：図8.4の点p_c，図8.4の点p_c，ニュートン法の初期点
double g_p_c[N];

// 関数f(x)の第1成分（double用）-----------------------------------------------------
double f_1_double(double x_1, double x_2) {
	return  2.0 * x_1 * x_1 * x_1 + 2. * x_1 * x_2 - 22.0 * x_1 + x_2 * x_2 + 13.0;
}
// 関数f(x)の第2成分（double用）-----------------------------------------------------
double f_2_double(double x_1, double x_2) {
	return  x_1 * x_1 + 2. * x_1 * x_2 + 2. * x_2 * x_2 * x_2 - 14. * x_2 + 9.;
}

// 関数f(x)の第1成分（自動微分型用）-------------------------------------------------
BU_N f_1(BU_N x_1, BU_N x_2) {
	return  2.0 * x_1 * x_1 * x_1 + 2. * x_1 * x_2 - 22.0 * x_1 + x_2 * x_2 + 13.0;
}
// 関数f(x)の第2成分（自動微分型用）-------------------------------------------------
BU_N f_2(BU_N x_1, BU_N x_2) {
	return  x_1 * x_1 + 2. * x_1 * x_2 + 2. * x_2 * x_2 * x_2 - 14. * x_2 + 9.;
}

// ホモトピー関数h(x)----------------------------------------------------------------
void func(BU_N h[], BU_N x[]) {

	BU_N BU_f[N - 1];

	// x[0]は式中のx_1,  x[1]：x_2 , x[2]： t

	// ホモトピー方程式
	h[0] = f_1(x[0], x[1]) + (x[2] - 1.0) * g_f_gamma_1;
	h[1] = f_2(x[0], x[1]) + (x[2] - 1.0) * g_f_gamma_2;

	// 平面q(x)=0の方程式
	h[2] = g_n[0] * (x[0] - g_p_c[0])
		+ g_n[1] * (x[1] - g_p_c[1])
		+ g_n[2] * (x[2] - g_p_c[2]);
}

// ニュートン法（多変数版）----------------------------------------------------------
void newton_n(void (*func)(BU_N [], BU_N []), double x[], int* r) {

	BU_N BU_x[N], BU_y[N];
	double f[N], J[N][N], delta_x[N];

	int iter = 0;
	do {
		// ベクトルxの値を自動微分型変数へ代入
		init(BU_x, x);

		// 自動微分型で関数値と偏微分値を計算
		func(BU_y, BU_x);

		// 関数値と偏微分値を分離して抽出
		split(BU_y, f, J);

		// 近似解の修正量△xを計算
		piboted_gauss(J, f, delta_x);

		// 近似解xの更新
		for (int i = 0; i < N; i++) x[i] -= delta_x[i];

		iter++;
	} while (vector_norm_max(delta_x) > EPS && iter < MAX_ITER);

	// 戻り値をセット
	*r = iter;
}

// ホモトピー法----------------------------------------------------------------------
int homotopy(void (*func)(BU_N [], BU_N []), double x_s[]) {

	// f(γ)の値をホモトピー関数h(x)にセットする
	g_f_gamma_1 = f_1_double(x_s[0], x_s[1]);
	g_f_gamma_2 = f_2_double(x_s[0], x_s[1]);

	// 解曲線上の追跡点（更新後）
	double p_new[N] = { x_s[0], x_s[1] , 0 };
	// 解曲線上の追跡点（更新前）
	double p_old[N] = { x_s[0], x_s[1], -DELTA_S };

	int i, iter_n, iter_h = 0, num = 0;

	while (fabs( p_new[2] ) < 5.0 && iter_h < MAX_STEP ) {

		// 解曲線の接線方向ベクトルを計算
		for (i = 0; i < N; i++) {
			g_n[i] = p_new[i] - p_old[i]; 
		}

		// Newton法の初期点p_cを計算
		for (i = 0; i < N; i++) {
			g_p_c[i] = p_new[i] + DELTA_S * g_n[i] / vector_norm_max(g_n);
		}

		// 次の反復の準備として更新後の追跡点を更新前の追跡点にセット
		for (i = 0; i < N; i++) p_old[i] = p_new[i];

		// p_newにNewton法の初期値をセット
		for (i = 0; i < N; i++) p_new[i] = g_p_c[i];

		// Newton法を実行 */
		newton_n(func, p_new, &iter_n);

		// ニュートン法の発散をチェック
		if (iter_n == MAX_ITER) {
			printf("解曲線の追跡に失敗しました．\n");
			return 1;
		}

		// 解曲線が平面t=0と交差していれば近似解を出力
		if ((p_new[2] - 1.0) * (p_old[2] - 1.0) < 0 || p_new[2] == 1.0) {
			printf("解 α_%d：\n", ++num );
			for (i = 0; i < N -1; i++) {
				printf("x[%d]=%lf\n", i, ( fabs(p_new[2] -1.0) *p_old[i] 
					+ fabs(1.0 -p_old[2]) *p_new[i] ) / fabs(p_new[2] -p_old[2]) );
			}
			printf("\n");
		}

		iter_h++;
	}

	if (iter_h == MAX_STEP) {
		printf("反復回数が上限に達しました．\n");
	} else {
		printf("tの値が設定範囲を越えました．\n");
	}

	return 0;
}

// メイン関数------------------------------------------------------------------------
int main() {

	double x_s[N - 1] = { -6.0, -6.0 };

	homotopy(func, x_s);

	return 0;
}