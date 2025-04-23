#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N   4           /* 次元数 */
#define EPS  1e-8       /* epsilonの設定 */
#define MAX_ITER 10     /* 最大反復回数 */

// 内積計算-------------------------------------------------------------------
double inner_product(const double x[], const double y[]) {

    double s = 0.0;
    for (int i = 0; i < N; i++) {
        s += x[i] * y[i];
    }

    return s;
}

// 最大値ノルム---------------------------------------------------------------
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

// 行列とベクトルの積---------------------------------------------------------
void matrix_vector_product(double c[], double a[][N], double b[]) {

    double s;
    int i, j;

    for (i = 0; i < N; i++){

        s = 0.0;
        for (j = 0; j < N; j++){

            s += a[i][j] * b[j];
        }
        c[i] = s;
    }
}

// 共役勾配法-----------------------------------------------------------------
void cg(double a[][N], double b[], double x[]) {

    double r[N], p[N], alpha, beta, tmp1[N], tmp2;
    int i, iter = 1;

    // x^(1) の計算------------------------------------------------------
    matrix_vector_product(tmp1, a, x);
    for (i = 0; i < N; i++) {
        r[i] = b[i] - tmp1[i];
        p[i] = r[i];
    }

    // A * p_k の計算 
    matrix_vector_product(tmp1, a, p);
    // (p_k, A*p_k) の計算 
    tmp2 = inner_product(p, tmp1);

    alpha = inner_product(p, p) / tmp2;

    for (i = 0; i < N; i++) x[i] = x[i] + alpha * p[i];
    for (i = 0; i < N; i++) r[i] = r[i] - alpha * tmp1[i];

    // x^(2) 以降の反復計算----------------------------------------------
    while (iter++ < MAX_ITER){

        beta = -inner_product(r, tmp1) / tmp2;
        for (i = 0; i < N; i++) p[i] = r[i] + beta * p[i];

        // A * p_k の計算 
        matrix_vector_product(tmp1, a, p);
        // (p_k, A*p_k) の計算 
        tmp2 = inner_product(p, tmp1);

        alpha = inner_product(p, r) / tmp2;

        for (i = 0; i < N; i++) x[i] = x[i] + alpha * p[i];
        for (i = 0; i < N; i++) r[i] = r[i] - alpha * tmp1[i];

        // 収束判定--------------------------------------------------
        if (vector_norm_max(r) <= EPS) {
            break;
        }
    }

    // 結果表示----------------------------------------------------------
    if (iter > MAX_ITER) {
        printf("収束しませんでした\n");
    }
    else {
        printf("解は以下の通り\n");
        printf("反復回数：%d回\n", iter);
        for (i = 0; i < N; i++) {
            printf("x_k[%d] = %lf\n", i, x[i]);
        }
    }
}

// メイン関数-----------------------------------------------------------------
int main(void){

    double a[N][N] = {
        { 6.0, 2.0, 1.0, 1.0 },
        { 2.0, 7.0, 3.0, 3.0 },
        { 1.0, 3.0, 8.0, 2.0 },
        { 1.0, 3.0, 2.0, 7.0 },
    };
    double b[N] = { 7.0, 6.0, -5.0, -3.0 };
    double x[N] = { 1.0, 1.0, 1.0, 1.0 };

    cg( a, b, x );

    return 0;
}



