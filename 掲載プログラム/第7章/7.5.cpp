#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#define EPS 1e-8       // 収束判定用閾値
#define MAX_ITER 7   // 最大反復回数

// 複素数クラス----------------------------------------------------------------
class Complex {

	// データメンバ------------------------------------------------------
	double Re, Im;

public:
	/* コンストラクタ */
	Complex(double x = 0, double y = 0) : Re(x), Im(y) { }

	/* 出力関数 */
	void print() {
		printf("%.6lf + %.6lfi ", Re, Im);
	}

	/* データへのアクセス */
	double get_Re() const { return Re; }
	double get_Im() const { return Im; }

	/* 絶対値の計算 */
	double abs() const {
		return sqrt(Re * Re + Im * Im);
	}

	// 演算子オーバーロード---------------------------------------------
	/* マイナス演算子 */
	Complex operator-() const {
		return Complex(-Re, -Im);
	}
	/* 加算 */
	friend Complex operator+(const Complex& x, const Complex& y) {
		return Complex(x.Re + y.Re, x.Im + y.Im);
	}
	/* 減算 */
	friend Complex operator-(const Complex& x, const Complex& y) {
		return Complex(x.Re - y.Re, x.Im - y.Im);
	}
	/* 乗算 */
	friend Complex operator*(const Complex& x, const Complex& y) {
		return Complex(x.Re * y.Re - x.Im * y.Im,
			x.Re * y.Im + x.Im * y.Re);
	}
	/* 除算 */
	friend Complex operator/(const Complex& x, const Complex& y) {
		double work = y.Re * y.Re + y.Im * y.Im;

		return Complex((x.Re * y.Re + x.Im * y.Im) / work,
			(-x.Re * y.Im + x.Im * y.Re) / work);
	}
	/* 指数関数 */
	friend Complex exp(const Complex& x) {
		double eR = std::exp(x.Re);
		return Complex(eR * cos(x.Im), eR * sin(x.Im));
	}
};

// 自動微分クラス（複素数版）---------------------------------------------
class BU_C {

	// データメンバ--------------------------------------------------
	Complex v;
	Complex d;

public:

	// コンストラクタ------------------------------------------------
	BU_C(double v_real = 0.0, double v_image = 0.0,
		double d_real = 0.0, double d_image = 0.0) {
		v = Complex(v_real, v_image);
		d = Complex(d_real, d_image);
	}

	BU_C(Complex x, Complex y = 0.0) {
		v = x;
		d = y;
	}

	/* 関数値，微分値の取得 */
	Complex value() const { return v; }
	Complex derivative() const { return d; }

	/* 微分値の初期化 */
	void init() {
		d = Complex(1.0, 0.0);
	}

	/* 出力関数 */
	void print() {
		v.print(); printf("< "); d.print(); printf(">");
	}

	// 演算子オーバーロード------------------------------------------
	/* 単項マイナス演算子 */
	BU_C operator-() const { return BU_C(-v, -d); }

	/* 加算 */
	friend BU_C operator+(const BU_C& x, const BU_C& y) {
		return BU_C(x.v + y.v, x.d + y.d);
	}
	/* 減算 */
	friend BU_C operator-(const BU_C& x, const BU_C& y) {
		return BU_C(x.v - y.v, x.d - y.d);
	}
	/* 乗算 */
	friend BU_C operator*(const BU_C& x, const BU_C& y) {
		return BU_C(x.v * y.v, x.d * y.v + x.v * y.d);
	}
	/* 除算 */
	friend BU_C operator/(const BU_C& x, const BU_C& y) {
		return BU_C(x.v / y.v,
			(x.d * y.v - x.v * y.d) / (y.v * y.v));
	}

	// 初等関数------------------------------------------------------
	friend BU_C exp(const BU_C& x) {
		Complex e = exp(x.v);
		return BU_C(e, e * x.d);
	}
};

// 関数の定義--------------------------------------------------------------
BU_C func(BU_C x) {

	return 3.0*x*x*x + x*x + 2.0;
}

// ニュートン法-----------------------------------------------------------------
Complex newton_C(BU_C(*f)(BU_C), Complex x_0, int* result, int* iterations) {

	BU_C x(x_0), y;
	Complex delta_x;
	int iter = 0, converged = 0;

	while (iter < MAX_ITER) {
		x.init();

		y = f(x);
		delta_x = y.value() / y.derivative();
		x = x - delta_x;

		iter++;

		// 収束判定--------------------------------------------------
		if (delta_x.abs() <= EPS) {
			converged = 1;
			break;
		}
	}

	// 戻り値のセット---------------------------------------------------
	*iterations = iter;
	*result = converged;
	return x.value();
}

// メイン関数-------------------------------------------------------------------
int main() {

	Complex y;
	int iter, result;

	y = newton_C(func, Complex(1.0, 1.0), &result, &iter);

	// 結果表示--------------------------------------------------------
	if (result == 1) {
		printf("反復回数は%d回です．\n", iter);
		printf("解：x = "); y.print(); printf("\n");
	} else {
		printf("収束しませんでした．\n");
	}

	return 0;
}