#include <stdio.h>

#define N 100   // 分割数

// 関数の定義--------------------------------------------------------------
double p(double x) {
    return -x;
}
double q(double x) {
    return 3.0;
}
double r(double x) {
    return  3.0 * x * x + 2 * x - 6;

}

// 差分法------------------------------------------------------------------
void fdm(double start, double end, double y_s, double y_e,
    int n, double x[], double y[]) {

    double alpha[N + 1], beta[N + 1], gamma[N + 1], b[N + 1];
    double m, h = (end - start) / N;
    int i;

    // 分点 x_i の計算----------------------------------
    for (i = 0; i <= N; i++) {
        x[i] = start + i * h;
    }

    // 係数行列の初期化---------------------------------
    for (i = 1; i <= N - 1; i++) {
        alpha[i] = 1 + h / 2.0 * p(x[i]);
        beta[i] = -2.0 - q(x[i]) * h * h;
        gamma[i] = 1.0 - h / 2.0 * p(x[i]);
        b[i] = h * h * r(x[i]);
    }
    b[1] -= alpha[1] * y_s;
    b[N - 1] -= gamma[N - 1] * y_e;

    // yの境界値を代入-----------------------------------
    y[0] = y_s;
    y[N] = y_e;

    // 前進消去------------------------------------------
    for (int i = 1; i <= N - 2; i++) {
        m = alpha[i + 1] / beta[i];
        beta[i + 1] = beta[i + 1] - m * gamma[i];
        b[i + 1] = b[i + 1] - m * b[i];
    }

    // 後退代入------------------------------------------
    y[N - 1] = b[N - 1] / beta[N - 1];
    for (int i = N - 2; i >= 1; i--) {
        y[i] = (b[i] - gamma[i] * y[i + 1]) / beta[i];
    }
}

// メイン関数--------------------------------------------------------------
int main() {
    double y[N + 1];
    double x[N + 1];
    int i;

    fdm(0.0, 3.0, 0.0, 6.0, N, x, y);

    // 結果を表示する------------------------------------
    printf("  x　        y      \n");
    printf("--------------------\n");
    for (i = 0; i <= N; i++) {
        printf("%5.2lf    %8.5lf\n", x[i], y[i]);
    }

    return 0;
}