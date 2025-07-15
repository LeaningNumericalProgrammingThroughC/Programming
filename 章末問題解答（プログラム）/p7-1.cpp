#include <stdio.h>
#include <math.h>

#define EPS 1e-8
#define MAX_ITER 30

// 関数の定義-------------------------------------------------------------------
double func(double x) {

	return x * x * sin(x) - 3.0 * x -2.0;
}

// 割線法-----------------------------------------------------------------------
double secant_method(double x_0, double x_1, int* r) {

	double y_0, y_1;
	double delta_x;
	int iter = 0;

	y_0 = func(x_0);

	do {

		y_1 = func(x_1);

		delta_x = y_1 * (x_1 - x_0) / (y_1 - y_0);

		x_0 = x_1;
		y_0 = y_1;
		x_1 = x_1 - delta_x;

		iter++;

	} while (fabs(delta_x) > EPS and iter < MAX_ITER);

	*r = iter;
	return x_1;
}

// メイン関数-------------------------------------------------------------------
int main() {

	double y;
	int iter;

	y = secant_method(1.0, 2.0, &iter);

	if (iter == MAX_ITER) {
		printf("収束しませんでした．\n");
	}
	else {
		printf("反復回数は%d回，解は%lfです．\n", iter, y);
	}

	return 0;
}


