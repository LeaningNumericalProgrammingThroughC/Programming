#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#define N 10
#define EPS 10e-8

// T_n(x)の値を計算
double cheby(double x, int n) {

    if (n == 0) return 1;
    if (n == 1) return x;
    return  2.0 * x * cheby(x, n - 1) - cheby(x, n - 2);
}

// チェビシェフ多項式の値を計算
double poly_cheby(double b[], double x) {

    double s;
    int i;

    s = 0.0;
    for (i = 0; i <= N; i++) s += b[i] * cheby(x, i);

    return s;
}

// シンプソン公式
double simpson_2(double a, double b, int n, double (*f)(double, int), int m) {

    double s, h;
    int i;

    h = (b - a) / (2.0 * n);

    s = f(a, m) + 4.0 * f(a + (2 * n - 1) * h, m) + f(b, m);
    for (i = 1; i < n; i++) {

        s += 4.0 * f(a + (2 * i - 1) * h, m) + 2.0 * f(a + 2 * i * h, m);
    }
    s *= h / 3.0;

    return s;
}

// f(cosθ)・T_m(θ) を計算する関数
double func_p11_3(double theta, int m) {

    double x = cos(theta);

    return (x * x * sin(3.0 * x) + 1.0) * cos((double)m * theta);
}

// 関数f(x)の定義
double func_p11_3_org(double x) {
    return x * x * sin(3.0 * x) + 1.0;
}

// メイン関数
int main() {

    double a = -1.0, b = 1.0, c[N + 1], x;
    int i, j, n = 50;
    double h = (b - a) / n;

    for (i = 0; i <= N; i++) {

        c[i] = simpson_2(0.0, M_PI, n, func_p11_3, i);

        if (i == 0) {
            c[i] = c[i] / M_PI;
        }
        else {
            c[i] = 2.0 / M_PI * c[i];
        }
    }

    // 結果表示
    printf("チェビシェフ多項式近似 ： %10.7f T_%d", c[N], N);
    for (i = N - 1; i >= 0; i--) {
        printf(" %c %10.7f T_%d", c[i] >= 0 ? '+' : '-', fabs(c[i]), i);
    }
    printf("\n\n");

    printf("       x              f(x)            P(x)\n");
    x = a;
    for (i = 0; i <= n; i++) {

        printf("%12.7f ,  %12.7f,   %12.7f\n", x, poly_cheby(c, x), func_p11_3_org(x));
        x += h;
    }

    return 0;
}