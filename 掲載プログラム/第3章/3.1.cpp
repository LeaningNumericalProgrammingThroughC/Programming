#include <stdio.h>

#define N 3   // 次元数

// ガウスの消去法---------------------------------------------
void gauss(double a[][N], double b[], double x[])
{
    double m, s;
    int i, j, k;

    // 前進消去----------------------------------------
    for (i = 0; i < N - 1; i++) {

        for (j = i + 1; j < N; j++) {

            m = a[j][i] / a[i][i];

            for (k = i + 1; k < N; k++) {
                a[j][k] = a[j][k] - m * a[i][k];
            }
            b[j] = b[j] - m * b[i];
        }
    }

    // 後退代入----------------------------------------
    x[N - 1] = b[N - 1] / a[N - 1][N - 1];
    for (i = N - 2; i >= 0; i--) {

        s = 0;
        for (j = i + 1; j < N; j++) {
            s += a[i][j] * x[j];
        }

        x[i] = (b[i] - s) / a[i][i];
    }
}

// メイン関数-------------------------------------------------
int main(){

    double a[N][N] = {
        {  2.0,  4.0, -2.0},
        {  3.0, -2.0,  1.0},
        { -2.0, -2.0,  3.0}
    };
    double b[N] = { 8.0, 8.0, -1.0 }, x[N];
    int i;

    gauss(a, b, x);

    printf("解は次の通りです．\n");
    for (i = 0; i < N; i++) {
        printf("x[%d] = %10.6lf\n", i, x[i]);
    }

    return 0;
}
