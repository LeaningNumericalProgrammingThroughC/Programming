#include <stdio.h>
#include <math.h>

#define N 100         // 分割数
#define EPS 1e-10     // 収束判定閾値
#define MAX_ITER 100   // 最大反復数

// 最大値ノルム------------------------------------------------------------
double vector_norm_max(double x[]) {

    double max;
    int i;

    max = fabs(x[0]);
    for (i = 1; i < N; i++) {
        if (fabs(x[i]) > max) {
            max = fabs(x[i]);
        }
    }

    return max;
}

// ニュートン法の修正量の計算----------------------------------------------
void delta_y(double alpha[], double beta[], double gamma[], double b[], double y[]) {

    double m;

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

    y[0] = y[N] = 0.0;
}

// 差分法------------------------------------------------------------------
void fdm(double start, double end, double y_s, double y_e, double y[], int* r) {

    double alpha[N + 1], beta[N + 1], gamma[N + 1], x, b[N + 1];
    double h = (end - start) / N, h_y = (y_e - y_s) / N;
    double d_y[N + 1], h_h, y_y, denominator;
    int i, iter = 0;

    // 初期近似解の計算---------------------------------
    for (i = 0; i <= N; i++) {
        y[i] = y_s + i * h_y;
    }

    // ニュートン法の反復計算---------------------------
    while (iter <= MAX_ITER) {

        // yの境界値を代入-----------------------------------
        y[0] = y_s;
        y[N] = y_e;

        // 差分方程式の各要素を計算--------------------------
        h_h = h * h;
        for (i = 1; i <= N - 1; i++) {

            x = start + i * h;

            y_y = y[i] * y[i];
            denominator = (x + 2.0) * (x + 2.0) * h;

            alpha[i] = 1.0 / h_h - y_y / denominator;
            beta[i] = -2.0 / h_h + 2.0 * y[i] * (y[i + 1] - y[i - 1]) / denominator;
            gamma[i] = 1.0 / h_h + y_y / denominator;
            b[i] = (y[i - 1] - 2.0 * y[i] + y[i + 1]) / h_h + y_y * (y[i + 1] - y[i - 1]) / denominator;
        }

        for (i = 1; i <= N - 1; i++) b[i] = -b[i];

        // ニュートン法の修正量を計算------------------------
        delta_y(alpha, beta, gamma, b, d_y);

        for (i = 1; i <= N - 1; i++) y[i] += d_y[i];

        // 収束判定------------------------------------------
        if (vector_norm_max(d_y) < EPS) break;

        iter++;
    }

    // 戻り値の設定
    *r = iter;
}

// メイン関数--------------------------------------------------------------
int main() {

    double y[N + 1], a = 0.0, b = 5.0;
    double y_a = 0.0, y_b = 3.0, h = (b - a) / N;
    int i, iter;

    // 差分法--------------------------------------------
    fdm(a, b, y_a, y_b, y, &iter);

    // 結果を表示する------------------------------------
    if (iter > MAX_ITER) {
        printf("ニュートン法が収束しませんでした．\n");
    }
    else {
        printf("解は以下の通り．（反復回数：%d回）\n", iter);
        printf("    x　        y      \n");
        printf("--------------------\n");
        for (i = 0; i <= N; i++) {
            printf("%7.2lf    %8.5lf\n", a+i*h, y[i]);
        }
    }

    return 0;
}