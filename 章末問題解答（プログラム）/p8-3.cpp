#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <cmath>

#define N  2            // 次元数
#define EPS 1e-8        // 収束判定用
#define MAX_ITER 100    // 最大反復回数
#define N_ROW 	13 * 13 // 収束状況を調べるための表の行数（初期点の個数）

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

// 関数の定義----------------------------------------------------------------------
void func(BU_N f[], BU_N x[]) {
	f[0] = 2.0 * x[0] * x[0] * x[0] + 2.0 * x[0] * x[1] - 22.0 * x[0] + x[1] * x[1] + 13.0;
	f[1] = x[0] * x[0] + 2.0 * x[0] * x[1] + 2.0 * x[1] * x[1] * x[1] - 14.0 * x[1] + 9.0;
}

// ニュートン法（多変数版）--------------------------------------------------------
void newton_n(void (*func)(BU_N f[], BU_N x[]), double x[], int* r) {

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

// メイン関数----------------------------------------------------------------------
int main() {

	const int n_row = N_ROW, n_col = 2 * N + 1;
	double x[N] = { 1.0,1.0 };
	
	// 見つかった解のリスト：solutions[i][j]は第i番目の解の第j成分
	double solutions[N_ROW][2];

	// ニュートンを反復実行させた各回の結果を保存するためのテーブル
	// 1列目:初期値x[0]，2列目:初期値x[1]，3列目:収束先αのα[0]，4列目:収束先αのα[1]
	// 5列目:収束先の解番号（第i番目の解に収束した場合はこの列の値はi）
	double table[n_row][n_col] = { 0.0 };
	
	int i, j, iter, row, num_sol;
	
	// 初期点を変更してニュートン表を実行．結果はtableに記録----------------
	row = 0;
	for (i = 0; i <= 12; i++) {
		for (j = 0; j <= 12; j++) {

			// 初期値をセット
			x[0] = -6.0 + (double)i;
			x[1] = -6.0 + (double)j;

			// 初期値を表に記入
			table[row][0] = x[0];
			table[row][1] = x[1];

			// ニュートン法を実行
			newton_n(func, x, &iter);

			// 収束値を表に記入
			if (iter == MAX_ITER) {
				table[row][N] = table[row][N+1] = NAN;
			} else {
				table[row][N] = x[0];
				table[row][N+1] = x[1];
			}

			row++;
		}
	}

	// 解のリストを作成しつつ，何番目の解に収束したかを表に記入------------
	num_sol = 0;
	for (i = 0; i < n_row; i++) {

		for (j = 0; j < num_sol; j++) {
			if (abs(table[i][N] - solutions[j][0]) < EPS
				&& abs(table[i][N+1] - solutions[j][1]) < EPS) {
				break;
			}
		}

		if (j == num_sol) {
			solutions[j][0] = table[i][N];
			solutions[j][1] = table[i][N+1];
			num_sol++;
		}

		table[i][2*N] = j;
	}

	// 表をcsvファイルに出力-----------------------------------------------
	FILE* fp = fopen("output.csv", "w");
	if (fp == NULL) {
		printf("ファイルを開けませんでした。\n");
		return 1;
	}

	//見出し行を出力
	for (i = 0; i < num_sol; i++) {
		fprintf(fp, ",");
		fprintf(fp, "(%lf，%lf)", solutions[i][0], solutions[i][1]);
	}
	fprintf(fp, "\n");


	// 表の各行を出力
	for (i = 0; i < n_row; i++) {

		fprintf(fp, "%5.2lf", table[i][0]);

		for (j = 0; j < num_sol; j++) {

			fprintf(fp, ",");
			if (j == table[i][2*N]) {
				fprintf(fp, "%5.2lf", table[i][1]);
			}
		}

		fprintf(fp, "\n");
	}

	// 解もプロットするためにファイルに出力
	for (i = 0; i < num_sol; i++) {

		fprintf(fp, "%5.2lf", solutions[i][0]);

		for (j = 0; j < num_sol; j++) {

			fprintf(fp, ",");
			if (j == i) {
				fprintf(fp, "%5.2lf", solutions[i][1]);
			}
		}

		fprintf(fp, "\n");
	}

	fclose(fp);

	return 0;
}