#include <stdio.h>
#define _USE_MATH_DEFINES  // Visual Studio等でM_PIを使う場合必要
#include <math.h>

#define N 100    // 標本点数

// 関数の定義----------------------------------------------------------------
double func(double x) {
    return x + cos(x);
}

// ガウス・チェビシェフ積分公式----------------------------------------------
double gauss_chebyshev(double (*f)(double), int n) {

    double s = 0.0, theta;
    int i;

    for (i = 1; i <= n; i++) {
        theta = (2.0 * i - 1.0) * M_PI / (2.0 * n);
        s += f(cos(theta)) * sin(theta);
    }

    return M_PI / n * s;
}

// メイン関数----------------------------------------------------------------
int main() {

    printf("ガウス・チェビシェフ積分公式: %12.10lf\n", gauss_chebyshev(func, N));

    return 0;
}