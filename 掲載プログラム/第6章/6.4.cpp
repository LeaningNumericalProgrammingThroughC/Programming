#include <stdio.h>
#include <cmath>

// 自動微分クラス（ボトムアップ型）--------------------------------------------
class BU {

	// データメンバ---------------------------------------------------
    double v;
    double d;

public:
	/* コンストラクタ */
	BU(double a = 0.0, double b = 0.0) : v(a), d(b) { }

    /* 関数値，微分値の取得 */
	double value() const { return v; }
	double derivative() const { return d; }

	/* 微分値の初期化 */
	void init() {
		d = 1.0;
	}

	/* 出力関数， 値 < 微分値 > の形式で表示 */
	void print() {
		printf("%lf < %lf >", v, d);
	}

	// 演算子オーバーロード-------------------------------------------
	/* 単項マイナス演算子 */
	BU operator-() const { return BU(-v, -d); }

	/* 加算 */
	friend BU operator+(const BU& x, const BU& y) {
		return BU(x.v + y.v, x.d + y.d);
	}
	/* 減算 */
	friend BU operator-(const BU& x, const BU& y) {
		return BU(x.v - y.v, x.d - y.d);
	}
	/* 乗算 */
	friend BU operator*(const BU& x, const BU& y) {
		return BU(x.v * y.v, x.d * y.v + x.v * y.d);
	}
	/* 除算 */
	friend BU operator/(const BU& x, const BU& y) {
		return BU( x.v / y.v, (x.d*y.v  - x.v*y.d) / (y.v*y.v) );
	}

	// 初等関数-------------------------------------------------------
	friend BU sqrt(const BU& x){
		double s = std::sqrt(x.v);
		return BU(s, 0.5 * x.d / s);
	}
	friend BU sin(const BU& x) {
		return BU(std::sin(x.v), std::cos(x.v) * x.d);
	}
	friend BU cos(const BU& x) {
		return BU(std::cos(x.v), -std::sin(x.v) * x.d);
	}
	friend BU exp(const BU& x){
		double e = std::exp(x.v);
		return BU(e, e * x.d);
	}
	friend BU log(const BU& x){
		return BU(std::log(x.v), x.d / x.v);
	}
};

// メイン関数------------------------------------------------------------------
int main() {

	BU x(2.0,1.0), y;

	y =  sin(x * x * x - exp(x)) / (x * x + 2.0);

	printf("y = ");	y.print(); printf("\n");

	return 0;
}