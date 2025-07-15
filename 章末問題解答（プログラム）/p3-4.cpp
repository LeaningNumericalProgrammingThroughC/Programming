#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 4          // 次元数 
#define EPS 1e-10    // 正則性判定用閾値 

// 逆行列計算-------------------------------------------------------------------
void inverse_matrix(double a[][N], double b[][N]) {

    int i, j, k, max_i;
    double m, max, tmp;

    // 前進消去------------------------------------------------------------
    for (i = 0; i < N; i++) {

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

            for (j = 0; j < N; j++) {
                tmp = b[i][j]; b[i][j] = b[max_i][j]; b[max_i][j] = tmp;
            }
        }

        // 前進消去-------------------------------------------------
        for (j = 0; j < N; j++) {

            if (i == j) continue;

            m = a[j][i] / a[i][i];

            for (k = i + 1; k < N; k++) {
                a[j][k] = a[j][k] - m * a[i][k];
            }

            for (k = 0; k < N; k++) {
                b[j][k] = b[j][k] - m * b[i][k];
            }

        }
    }

    // 対角成分で割って解を求める------------------------------------------
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            b[i][j] = b[i][j] / a[i][i];
        }
    }
}

// メイン関数----------------------------------------------------------------------
int main(void) {

    double a[N][N] = {
        {  0.0, -3.0, -2.0, -4.0},
        {  3.0, -2.0,  2.0, -1.0},
        { -2.0,  0.0,  3.0,  2.0},
        {  5.0,  2.0, -3.0,  3.0}
    };

    double b[N][N] = {
        {  1.0,  0.0,  0.0,  0.0},
        {  0.0,  1.0,  0.0,  0.0},
        {  0.0,  0.0,  1.0,  0.0},
        {  0.0,  0.0,  0.0,  1.0}
    };

    int i, j;

    inverse_matrix(a, b);

    printf("逆行列：\n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%10.6lf ", b[i][j]);
        }
        printf("\n");
    }

    return 0;
}
