#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N    4       // 次元数 
#define EPS 1e-10    // 正則性判定用

// LU分解----------------------------------------------------------------------
void lu_decomposition(double a[][N], int p[]) {

    double m, s, max, tmp;
    int i, j, k, max_i;

    for (i = 0; i < N - 1; i++)
    {
        // 部分ピボット選択-----------------------------------
        // 最大行の検索
        max = fabs(a[i][i]);
        max_i = i;
        for (j = i + 1; j < N; j++) {

            if (fabs(a[j][i]) > max) {
                max = fabs(a[j][i]);
                max_i = j;
            }
        }

        // 入れ換えた行を配列に記憶
        p[i] = max_i;
        
        // 正則性の判定
        if (max < EPS) {
            printf("係数行列が正則ではありません");
            exit(1);
        }
        
        // 行の入れ換え
        if (max_i != i) {

            for (j = 0; j < N; j++) {
                tmp = a[i][j]; a[i][j] = a[max_i][j]; a[max_i][j] = tmp;
            }
        }
      
        // 前進消去-------------------------------------------
        for (j = i + 1; j < N; j++) {

            m = a[j][i] / a[i][i];

            // 行列Uの要素を計算
            for (k = i + 1; k < N; k++) {
                a[j][k] = a[j][k] - m * a[i][k];
            }

            // 行列Lの要素を格納
            a[j][i] = m;
        }
    }
}

// 前進代入と後退代入----------------------------------------------------------
void lu_substitution(double a[][N], double b[], int p[]) {

    double tmp, s;
    int i, j, k;

    // 方程式の右辺Pbを求める-----------------------------------
    for (i = 0; i < N-1; i++) {
        tmp = b[i]; b[i] = b[p[i]]; b[p[i]] = tmp;
    }

    // 前進代入-------------------------------------------------
    for (i = 0; i < N; i++) {

        for (j = 0; j < i; j++) {
            b[i] -= a[i][j] * b[j];
        }
    }

    // 後退代入-------------------------------------------------
    b[N - 1] = b[N - 1] / a[N - 1][N - 1];
    for (i = N - 2; i >= 0; i--) {

        s = 0;
        for (j = i + 1; j < N; j++) {
            s += a[i][j] * b[j];
        }
        b[i] = (b[i] - s) / a[i][i];
    }
}

// メイン関数------------------------------------------------------------------
int main(void){

    double a[N][N] = {
        {  2.0,  2.0,  1.0,  2.0 },
        {  2.0,  3.0, -1.0,  3.0 },
        { -3.0,  4.0,  0.0,  3.0 },
        {  1.0,  3.0, -2.0,  1.0 }
    };
    double b[N] = { 7.0, 6.0, -5.0, -3.0 };
    int i, p[N];

    lu_decomposition(a, p);

    lu_substitution(a, b, p);

    // 結果出力-------------------------------------------------
    printf("解は次の通り\n");
    for (i = 0; i < N; i++) {
        printf("x[%d] = %10.6lf\n", i, b[i]);
    }

    return 0;
}