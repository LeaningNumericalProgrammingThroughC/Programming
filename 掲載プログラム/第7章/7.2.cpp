#include <stdio.h>
#include <math.h>

#define EPS 1e-8         // 収束判定用閾値
#define MAX_ITER 50      // 最大反復回数

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
BU func(BU x) {

	return x * x * x - 3 * x * x - exp(2 * x - 5) + 3;
}

// ニュートン法-----------------------------------------------------------------
double newton( BU (*f)(BU), double x_0, int *result, int* iterations) {

	BU x(x_0), y;
	double delta_x;
	int iter = 0, converged = 0;

	while (iter < MAX_ITER) {

		// 微分値を初期化（0をセット）
		x.init();  

		y = f(x);
		delta_x = y.value() / y.derivative();
		x = x - delta_x;
		
		iter++;

		// 収束判定
		if (fabs(delta_x) <= EPS) {
			converged = 1; 
			break;
		}
	}

	// 戻り値のセット
	*iterations = iter;
	*result = converged;
	return x.value();
}

// メイン関数-------------------------------------------------------------------
int main() {

	double y;
	int iter, result;
	
	y = newton(func, 1.0, &result, &iter);
	
	// 結果表示--------------------------------------------------------
	if (result == 1) {
		printf("反復回数は%d回，解は%lfです．\n", iter, y);
	} else {
		printf("収束しませんでした．\n");
	}

	return 0;
}