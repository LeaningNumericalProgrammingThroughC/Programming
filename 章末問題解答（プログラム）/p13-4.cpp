#include <stdio.h>

// 被積分関数--------------------------------------------------------------------------
double f(double x, double y) {
	return x * x * y + x + 1.0;
}

// 二重積分----------------------------------------------------------------------------
double double_integral(double (*f)(double, double), double x_a, double x_b, int x_n,
	               double y_a, double y_b, int y_n) {

	int i, j, x_m = 2 * x_n, y_m = 2 * y_n;
	double x, y, x_s, y_s;
	double x_h = (x_b - x_a) / x_m, y_h = (y_b - y_a) / y_m, alpha;

	// 二重積分の外側の積分を計算する-----------------------------------
	x_s = 0.0;
	for (i = 0; i < x_m; i++) {

		x = x_a + i * x_h;

		// 二重積分の内側の積分を計算する-----------------------
		y_s = 0.0;
		for (j = 0; j < y_m; j++) {

			y = y_a + j * y_h;

			if (j == 0 || j == y_m) alpha = 1.0;
			else if (j % 2 == 0)  alpha = 2.0;
			else alpha = 4.0;

			y_s += alpha * f(x, y);
		}
		y_s *= y_h / 3.0;
		
		if (i == 0 || i == x_m) alpha = 1.0;
		else if (i % 2 == 0)  alpha = 2.0;
		else alpha = 4.0;

		x_s += alpha * y_s;
	}
	x_s *= x_h / 3.0;

	return x_s;
}

// メイン関数-------------------------------------------------------------------------
int main() {
	
	printf("積分値：%8.6lf\n", double_integral(f, 0.0, 1.0, 100, 0.0, 1.0, 100));

	return 0;
}