#include <stdio.h>
#include <math.h>

// 関数の定義-----------------------------------------------------------
double func(double x) {
    return x + cos(x);
}

// 台形公式-------------------------------------------------------------
double trapezoidal(double (*f)(double), double a, double b, int n) {

    double s, h = (b - a) / n;
    int i;

    s = f(a) + f(b);
    for (i = 1; i < n; i++) {
        s += 2.0 * f(a + i * h);
    }

    return s * h / 2.0;
}

// シンプソン公式-------------------------------------------------------
double simpson(double (*f)(double), double a, double b, int n) {

    double s, h = (b - a) / (2.0 * n);
    int i;

    s = f(a) + 4.0*f(a + (2*n - 1)*h) + f(b);
    for (i = 1; i < n; i++) {
        s += 4.0 * f(a + (2*i-1)*h) + 2.0 * f(a + 2*i*h);
    }

    return s * h / 3.0;
}

// メイン関数-----------------------------------------------------------
int main() {

    double a = -1, b = 1;
    int n = 100;

    // 台形公式
    printf("台形公式　　　 ： %12.10lf\n", trapezoidal(func, a, b, n));
    // シンプソン公式
    printf("シンプソン公式 ： %12.10lf\n", simpson(func, a, b, n));

    return 0;
}
