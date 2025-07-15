#include <stdio.h>
#include <math.h>

#define MAX_ITER 10       // 最大反復回数

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
void trapezoidal_romberg(double (*f)(double), double a, double b, int n, double exact_val) {

    double T[MAX_ITER][MAX_ITER], s, h = (b - a) / n;
    int i, j, n_org = n;

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

        // 次の反復における分点数を計算
        n = 2 * n;
    }

    // 結果出力---------------------------------------------------------------------------
    h = (b - a) / n_org;
    printf("     区間幅           台形公式の誤差      ロンバーグ積分の誤差\n");
    printf("---------------------------------------------------------------\n");
    for (i = 0; i < MAX_ITER; i++) {
        printf("h = %10.5e   %20.16lf   %20.16lf\n", pow(0.5,i)*h, 
                     fabs(T[i][0] - exact_val), fabs(T[i][i] - exact_val));
    }
}

// メイン関数-----------------------------------------------------------------------
int main() {

    double a = -1.0, b = 1.0;
    double exact_val = sin(1.0) - sin(-1.0);
    
    trapezoidal_romberg(func, a, b, 10, exact_val);

    return 0;
}
