#include <stdio.h>

#define N 100    // 分割数

// 微分方程式 dy/dx = f(x, y) を定義--------------------------------
double f(double x, double y) {

    return 8 * x * y;
}

// リープフロッグ法-------------------------------------------------------
void leapfrog(double a, double b, double y[], int n,
    double (*f) (double, double)) {

    double x = a, h = (b - a) / n;

    y[1] = y[0] + h * f(x, y[0]);

    x = x + h;
    for (int i = 1; i < n; i++) {
        y[i + 1] = y[i - 1] + 2.0 * h * f(x, y[i]);
        x = x + h;
    }
}

// メイン関数-------------------------------------------------------
int main() {

    double a = 0.0, b = 1.0, h = (b - a) / N, y[N + 1];
    int i;

    // リープ・フロッグ法
    y[0] = 1.0;
    leapfrog(a, b, y, N, f);
    
    // 結果出力---------------------------------------------
    printf("  x　        y      \n");
    printf("--------------------\n");
    for (i = 0; i <= N; i++) {
        printf("%5.2lf    %8.5lf\n", i * h, y[i]);
    }

    return 0;
}
