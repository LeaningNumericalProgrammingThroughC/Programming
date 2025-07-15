#include <stdio.h>
#include <math.h>

#define EPS 1e-8
#define MAX_ITER 100

// 自動微分クラス（ボトムアップ型）--------------------------------------------
class BU {

	double v;
	double d;

public:

	/* コンストラクタ */
	BU(double a = 0.0, double b = 0.0) {
		v = a; d = b;
	}
	/* 関数値，微分値の取得 */
	double value() const { return v; }
	double derivative() const { return d; }

	/* 微分値の初期化 */
	void init() {
		d = 1.0;
	}
	/* 出力関数 */
	void print() {
		printf("%lf < %lf >", v, d);
	}

	/* マイナス演算子 */
	BU operator-() const { return BU(-v, -d); }

	/* 加算 (BU ＋ BU) */
	friend BU operator+(const BU& x, const BU& y) {
		return BU(x.v + y.v, x.d + y.d);
	}
	/* 減算 (BU - BU) */
	friend BU operator-(const BU& x, const BU& y) {
		return BU(x.v - y.v, x.d - y.d);
	}
	/* 乗算 (BU * BU) */
	friend BU operator*(const BU& x, const BU& y) {
		return BU(x.v * y.v, x.d * y.v + x.v * y.d);
	}
	/* 除算 (BU / BU) */
	friend BU operator/(const BU& x, const BU& y) {
		return BU(x.v / y.v, (x.d * y.v - x.v * y.d) / (y.v * y.v));
	}

	/* 初等関数 */
	friend BU sqrt(const BU& x) {
		return BU(sqrt(x.v), 0.5 * x.d / (sqrt(x.v)));
	}
	friend BU sin(const BU& x) {
		return BU(sin(x.v), cos(x.v) * x.d);
	}
	friend BU cos(const BU& x) {
		return BU(cos(x.v), -sin(x.v) * x.d);
	}
	friend BU exp(const BU& x) {
		return BU(exp(x.v), exp(x.v) * x.d);
	}
	friend BU log(const BU& x) {
		return BU(log(x.v), x.d / x.v);
	}

};

// 関数の定義-------------------------------------------------------------------
BU func_BU(BU x) {

	return x * x * sin(x) - 3.0 * x - 2.0;
}

// 関数の定義-------------------------------------------------------------------
double func_double(double x) {

	return x * x * sin(x) - 3.0 * x - 2.0;
}

// 簡易ニュートン法-------------------------------------------------------------
double newton(double x, int* r) {

	BU x_BU(x), y;
	double df, delta_x;
	int iter = 0;

	x_BU.init();  // 微分値を初期化（1をセット）
	y = func_BU(x_BU);
	df = y.derivative();

	do {

		delta_x = func_double(x) / df;
		x = x - delta_x;

		iter++;

	} while (fabs(delta_x) > EPS and iter < MAX_ITER);

	*r = iter;
	return x;
}

// メイン関数-------------------------------------------------------------------
int main() {

	double y;
	int iter;

	y = newton(0.0, &iter);

	if (iter == MAX_ITER) {
		printf("収束しませんでした．\n");
	}
	else {
		printf("反復回数は%d回，解は%lfです．\n", iter, y);
	}

	return 0;
}