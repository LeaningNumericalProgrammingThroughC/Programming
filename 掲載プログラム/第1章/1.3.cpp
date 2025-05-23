#include <stdio.h>
#include <math.h>

int main() {
	float a, b, c, D, r;
	double da, db, dc, dD;

	// 係数の初期化
	a = 1.26356; b = 17834.6; c = 2.51522;

	// 判別式の計算
	D = b * b - 4.0 * a * c;


	printf("通常の解の公式を用いて計算\n");
	printf("x_1 = %.10f\n", (-b + sqrt(D)) / (2.0 * a));
	printf("x_2 = %.10f\n\n", (-b - sqrt(D)) / (2.0 * a));


	printf("有理化した公式を用いて計算\n");
	if (b >= 0) {
		printf("x_1 = %.10f\n", -2.0 * c / (b + sqrt(D)));
		printf("x_2 = %.10f\n\n", (-b - sqrt(D)) / (2.0 * a));
	}
	else {
		printf("x_1 = %.10f\n\n", (-b + sqrt(D)) / (2.0 * a));
		printf("x_2 = %.10f\n", -2.0 * c / (b - sqrt(D)));
	}

	da = (double)a; db = (double)b; dc = (double)c;
	dD = db * db - 4.0 * da * dc;
	printf("倍精度で計算\n");
	printf("x_1 = %.10f\n", (-db + sqrt(dD)) / (2.0 * da));
	printf("x_2 = %.10f\n\n", (-db - sqrt(dD)) / (2.0 * da));
}

