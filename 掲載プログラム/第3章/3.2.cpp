#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 4          // 次元数 
#define EPS 1e-10    // 正則性判定用閾値 

// 部分ピボット選択付きガウスの消去法----------------------------------------------
void pivoted_gauss(double a[][N], double b[], double x[])
{
    int i, j, k, max_i;
    double m, s, max, tmp;

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

// メイン関数----------------------------------------------------------------------
int main(void){

    double a[N][N] = {
        {  0.0, -3.0, -2.0, -4.0 },
        {  3.0, -2.0,  2.0, -1.0 },
        { -2.0,  0.0,  3.0,  2.0 },
        {  5.0,  2.0, -3.0,  3.0 }
    };
    double b[N] = { 2.0, -3.0, -2.0, 1.0 }, x[N];
    int i;

    pivoted_gauss(a, b, x);

    printf("解は次の通り\n");
    for (i = 0; i < N; i++) {
        printf("x[%d] = %10.6lf\n", i, x[i]);
    }

    return 0;
}

