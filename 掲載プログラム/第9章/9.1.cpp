#include <stdio.h>
#include <math.h>

#define N    4          // 次元数
#define EPS    1e-10    // 収束判定用
#define MAX_ITER  500   // 最大反復回数

// 行列とベクトルの乗算-----------------------------------------------
void mat_vec_product(double r[], double a[][N], double x[]) {
    int i, j;

    for (i = 0; i < N; i++) {
        r[i] = 0.0;
        for (j = 0; j < N; j++) {
            r[i] += a[i][j] * x[j];
        }
    }
}

// ベクトルの内積-----------------------------------------------------
double inner_product(double a[N], double b[N]) {
    int i;
    double s = 0.0;

    for (i = 0; i < N; i++) {
        s += a[i] * b[i];
    }

    return s;
}

// べき乗法-----------------------------------------------------------
double power_method(double A[][N], double x[], int* r) {

    double Ax[N], lambda, x_sq, norm_x;
    int i, iter = 0;

    do {

        // 行列Aとベクトルxの積を計算する
        mat_vec_product(Ax, A, x);
        
        // レイリー商により近似固有値を求める
        lambda = inner_product(Ax, x);
        
        // Axの2ノルムを計算する 
        x_sq = inner_product(Ax, Ax);
        norm_x = sqrt(x_sq);

        // 次の反復のために，Axの値をxに戻す
        for (i = 0; i < N; i++) x[i] = Ax[i] / norm_x;

        iter++;

    } while (fabs(x_sq - lambda*lambda) >= EPS && iter < MAX_ITER);

    *r = iter;
    return lambda;
}

// メイン関数---------------------------------------------------------
int main() {

    double A[N][N] = {
        {  8.0,  3.0, -1.0, -3.0},
        {  2.0, -1.0,  3.0,  4.0},
        {  2.0, -6.0,  2.0,  1.0},
        { -1.0, -2.0,  7.0, -2.0}
    };
    double x[N] = { 1.0, 1.0, 1.0, 1.0 }, lambda;
    int i, iter;

    lambda = power_method(A, x, &iter);

    // 結果の表示------------------------------------------
    if (iter == MAX_ITER) {
        printf("べき乗法が収束しませんでした．");
    } else {
        printf("反復回数は%d回です．\n", iter);
        printf("最大固有値は %f です．\n", lambda);
        printf("固有ベクトルは以下です，\n");
        for (i = 0; i < N; i++) {
            printf("x[%d] = %f\n",i, x[i]);
        }
    }
}
