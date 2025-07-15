#include <stdio.h>
#include <math.h>

#define N 4   /* 次元数 */

// コレスキー分解---------------------------------------------------------
void cholesky_decomp(double a[][N]) {

    double s;
    int i, j, k;

    for (i = 0; i < N; i++) {

        // 対角要素を求める------------------------------------
        s = 0.0;
        for (k = 0; k < i; k++) {
            s += a[i][k] * a[i][k];
        }
        a[i][i] = sqrt(a[i][i] - s);

        // 第i列の対角要素より下部分を求める-------------------
        for (j = i + 1; j < N; j++) {
            s = 0.0;
            for (k = 0; k < i; k++) {
                s += a[i][k] * a[j][k];
            }
            a[j][i] = (a[j][i] - s) / a[i][i];
        }
    }
}

// 前進代入と後退代入---------------------------------------------------------
void cholesky_substitution(double a[][N], double b[N]) {

    double s;
    int i, j;

    // Ly = bの解yを求める（yの値はbに格納）------------------------
    for (i = 0; i < N; i++) {

        s = 0.0;
        for (j = 0; j < i; j++) {

            s += a[i][j] * b[j];
        }
        b[i] = (b[i] - s) / a[i][i];
    }

    // L^tx = y の解xを求める（xの値はbに格納)------------------------
    for (i = N - 1; i >= 0; i--) {

        s = 0.0;
        for (j = i + 1; j < N; j++) {

            s += a[j][i] * b[j];
        }
        b[i] = (b[i] - s) / a[i][i];
    }
}

// メイン関数-----------------------------------------------------------------
int main(void) {

    double a[N][N] = {
        { 15.0,  7.0, -6.0,  8.0},
        {  7.0,  6.0, -7.0,  3.0},
        { -6.0, -7.0, 14.0,  1.0},
        {  8.0,  3.0,  1.0, 13.0}
    };
    double b[N] = { -2.0, 6.0, -3.0, 1.0 };
    int i;

    cholesky_decomp(a);

    cholesky_substitution(a, b);

    printf("解は次の通り\n");
    for (i = 0; i < N; i++) {
        printf("x[%d] = %8.5lf\n", i, b[i]);
    }

    return 0;
}