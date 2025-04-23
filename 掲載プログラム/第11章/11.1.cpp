#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M 8          // データ点数
#define N 4          // 多項式の未知係数の数（次数n + 1）
#define EPS 1e-10    // 正則性判定用閾値 

// 部分ピボット選択付きガウスの消去法----------------------------------------------
void pivoted_gauss(double a[][N], double b[], double x[]){

    double m, s, max, tmp;
    int i, j, k, max_i;

    // 前進消去------------------------------------------------------------
    for (i = 0; i < N - 1; i++)
    {
        // 部分ピボット選択-----------------------------------------
        max = fabs(a[i][i]);
        max_i = i;
        for (j = i + 1; j < N; j++) {

            if (fabs(a[j][i]) > max) {
                max = fabs(a[j][i]);
                max_i = j;
            }
        }

        if (max < EPS) {

            printf("係数行列が正則ではありません");
            exit(1);
        }
        else if (max_i != i) {

            for (j = i; j < N; j++) {
                tmp = a[i][j]; a[i][j] = a[max_i][j]; a[max_i][j] = tmp;
            }
            tmp = b[i]; b[i] = b[max_i]; b[max_i] = tmp;
        }

        // 前進消去-------------------------------------------------
        for (j = i + 1; j < N; j++) {

            m = a[j][i] / a[i][i];

            for (k = i + 1; k < N; k++) {
                a[j][k] = a[j][k] - m * a[i][k];
            }
            b[j] = b[j] - m * b[i];
        }
    }

    // 後退代入------------------------------------------------------------
    x[N - 1] = b[N - 1] / a[N - 1][N - 1];
    for (i = N - 2; i >= 0; i--) {

        s = 0;
        for (j = i + 1; j < N; j++) {
            s += a[i][j] * x[j];
        }
        x[i] = (b[i] - s) / a[i][i];
    }
}

// 最小二乗近似--------------------------------------------------------
void least_squares(double x[], double y[], double c[]) {

    double a[N][N], b[N];
    int i, j, k;

    // 行列の各成分を計算----------------------------------------
    for (i = 0; i < N; i++) {
        for (j = 0; j <= i; j++) {
            a[i][j] = 0.0;
            for (k = 0; k < M; k++) {
                a[i][j] += pow(x[k], i + j);
            }
            a[j][i] = a[i][j];
        }
    }

    // 右辺ベクトルを計算----------------------------------------
    for (i = 0; i < N; i++) {
        b[i] = 0.0;
        for (k = 0; k < M; k++) {
            b[i] += y[k] * pow(x[k], i);
        }
    }

    // ガウスの消去法で多項式の係数を求める----------------------
    pivoted_gauss(a, b, c);
}

// メイン関数----------------------------------------------------------
int main(void) {

    double x[M] = { 0.12, 0.25, 0.31, 0.33, 0.51, 0.72, 0.89, 0.91 };
    double y[M] = { 0.43, 1.12, 3.62, 4.57, 4.81, 6.22, 5.11, 3.86 };
    double c[N];
    int i;
    
    least_squares(x, y, c);

    // 結果出力--------------------------------------------------
    printf("最小二乗近似多項式の係数:\n");
    for (i = 0; i < N; i++) {
        printf("c[%d] = %10.6lf\n", i, c[i]);
    }

    return 0;
}
