#include <stdio.h>
#include <math.h>

#define EPS 1e-8          // 収束判定用閾値
#define MAX_ITER 20       // 最大反復回数

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

// 関数の定義-------------------------------------------------------------------------
double func(double x) {
    return x + cos(x);
}

// ロンバーグ積分-------------------------------------------------------------------
double romberg(double (*f)(double), double a, double b, int n) {

    double T[MAX_ITER][MAX_ITER], s, h = (b - a) / n;
    int i, j;

    T[0][0] = trapezoidal(f, a, b, n);

    for (i = 1; i < MAX_ITER; i++) {

        // 刻み幅を半分した際の面積を台形公式で計算-----------------------------
        h /= 2.0;
        s = 0.0;
        for (int k = 1; k <= n; k++) {
            s += f(a + (2 * k - 1) * h);
        }
        T[i][0] = 0.5 * T[i - 1][0] + h * s;

        // T表のi行目を右に向かって計算する-------------------------------------
        for (j = 1; j <= i; j++) {
            T[i][j] = (pow(4.0, j) * T[i][j - 1] - T[i - 1][j - 1]) / (pow(4, j) - 1.0);
        }

        // 収束判定-------------------------------------------------------------
        if (fabs(T[i][i] - T[i - 1][i - 1]) < EPS) {
            return T[i][i];
        }

        // 次の反復における分点の増加分を計算
        n = 2 * n;
    }
    return T[MAX_ITER - 1][MAX_ITER - 1];
}

// メイン関数-----------------------------------------------------------------------
int main() {

    double a = -1.0, b = 1.0;

    printf("ロンバーグ積分による積分値: %20.18lf\n", romberg(func, a, b, 10));

    return 0;
}
