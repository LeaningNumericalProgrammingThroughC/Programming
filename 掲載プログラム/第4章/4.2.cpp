#include <stdio.h>

#define N 4   /* 次元数 */

// 修正コレスキー分解---------------------------------------------------------
void cholesky_decomp(double a[][N]) {

    double s;
    int i, j, k;

    for (i = 0; i < N; i++) {

        // 対角行列Dを求める (行列Aの対角成分に格納）--------------------
        s = 0.0;
        for (k = 0; k < i; k++) {
            s += a[i][k] * a[k][k] * a[i][k];
        }
        a[i][i] -= s;

        // 下三角行列Lを求める（行列Aの下三角部分に格納）----------------
        for (j = i + 1; j < N; j++){

            s = 0.0;
            for (k = 0; k < i; k++){

                s += a[i][k] * a[k][k] * a[j][k];
            }
            a[j][i] = (a[j][i] - s) / a[i][i];
        }
    }
}

// 前進代入と後退代入---------------------------------------------------------
void cholesky_substitution(double a[][N], double b[N]) {

    double s;
    int i, j;

    // LD y = bの解yを求める（yの値はbに格納）------------------------
    for (i = 0; i < N; i++){

        s = 0.0;
        for (j = 0; j < i; j++){

            s += a[j][j] * a[i][j] * b[j];
        }
        b[i] = (b[i] - s) / a[i][i];
    }

    // L^t x = y の解xを求める（xの値はbに格納)------------------------
    for (i = N - 2; i >= 0; i--) {

        s = 0.0;
        for (j = i + 1; j < N; j++){

            s += a[j][i] * b[j];
        }
        b[i] -= s;
    }
}

// メイン関数-----------------------------------------------------------------
int main(void) {

    double a[N][N] = {
        {  3.0, -2.0,  1.0,  2.0},
        { -2.0,  6.0,  0.0, -2.0},
        {  1.0,  0.0,  5.0,  4.0},
        {  2.0, -2.0,  4.0,  6.0}
    };
    double b[N] = { 7.0, -6.0, 9.0, 6.0 };
    int i;

    cholesky_decomp(a);

    cholesky_substitution(a, b);

    printf("解は次の通り\n");
    for (i = 0; i < N; i++) {
        printf("x[%d] = %lf\n", i, b[i]);
    }

    return 0;
}