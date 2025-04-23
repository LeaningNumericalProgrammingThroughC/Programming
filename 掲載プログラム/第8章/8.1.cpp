#include <stdio.h>
#include <cmath>

#define N 3    // 次元数

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
	for (int i = 0; i < N; i++) { r.d[i] = a.d[i]*b.v + a.v*b.d[i]; }
	return r;
}
/* 除算 */
BU_N operator/(const BU_N& a, const BU_N& b) {
	BU_N r(a.v / b.v);
	for (int i = 0; i < N; i++) { 
		r.d[i] = ( a.d[i]*b.v - a.v*b.d[i] ) / (b.v * b.v); 
	}
	return r;
}

// 初等関数--------------------------------------------------------------------
BU_N sin(const BU_N& x) {
	BU_N r(std::sin(x.v));
	for (int i = 0; i < N; i++) { r.d[i] = std::cos(x.v) * x.d[i]; }
	return r;
}

   // ・・・　（以下，cos(x),exp(x)等も同様）


// メイン関数------------------------------------------------------------------
int main() {

	double w[N] = { 1,0,0 };

	BU_N x(2.0),y(3.0),z(4.0);

	/* 各変数の偏微分値の初期化 */
	x.d[0] = y.d[1] = z.d[2] = 1.0;

	BU_N f =  5.0 * x * x - 2.0 * x * y * y + 3.0 * y * z;
	
	printf("f <df/dx, df/dy, df/dz> = ");	f.print(); printf("\n");

	return 0;
}
