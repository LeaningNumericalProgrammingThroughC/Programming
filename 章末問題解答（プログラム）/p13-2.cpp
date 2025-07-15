#include <stdio.h>
#include <math.h>

// 関数の定義-----------------------------------------------------------
double func(double x) {
    return x + cos(x);
}

// 中点公式-------------------------------------------------------------
double midpoint(double (*f)(double), double a, double b, int n) {

    double s, h = (b - a) / (2.0 * n);
    int i;

    s = 0.0;
    for (i = 1; i <= n; i++) {
        s += f(a + (2.0 * i - 1.0) * h);
    }

    return 2.0 * h * s;
}

// メイン関数-----------------------------------------------------------
int main() {

    double a = -1, b = 1;
    int n = 100;

    // 中点公式
    printf("中点公式による積分値：%12.10lf\n", midpoint(func, a, b, n));

    return 0;
}
