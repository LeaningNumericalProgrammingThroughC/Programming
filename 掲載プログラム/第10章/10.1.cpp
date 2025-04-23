#include <stdio.h>
#include <math.h>

#define   N   5     // 次元数

// 内積計算------------------------------------------------------------------
double inner_product(const double x[], const double y[]) {

    double s = 0.0;
    for (int i = 0; i < N; i++) {
        s += x[i] * y[i];
    }

    return s;
}

// 2ノルム-------------------------------------------------------------------
double vector_norm2(double x[]) {

    double s;
    int i;

    s = 0.0;
    for (i = 0; i < N; i++) {
        s += x[i] * x[i];
    }

    return sqrt(s);
}

// 行列とベクトルの積---------------------------------------------------------
void matrix_vector_product(double c[], double a[][N], double b[]) {

    double s;
    int i, j;

    for (i = 0; i < N; i++) {

        s = 0.0;
        for (j = 0; j < N; j++) {

            s += a[i][j] * b[j];
        }
        c[i] = s;
    }
}

// 符号関数（ただし，0の場合は1を返す）---------------------------------
int sign(double x) {
    return (x >= 0.0) ? 1 : -1;
}

// 鏡映変換による三重対角化---------------------------------------------
void householder_trans(double a[][N]) {

    double s, w, alpha, u[N], p[N], q[N];
    int i, j, k;

    // i: 反射変換を適用する列（0～N-3まで）
    for (i = 0; i < N - 2; i++) {

        // uの構成---------------------------------------------
        s = 0.0;
        for (j = i + 1; j < N; j++) {
            s += a[j][i] * a[j][i];
        }
        alpha = -sign(a[i + 1][i]) * sqrt(s);

        for (j = 0; j <= i; j++) u[j] = 0.0;
        u[i + 1] = a[i + 1][i] - alpha;
        for (j = i + 2; j < N; j++) {
            u[j] = a[j][i];
        }

        w = vector_norm2(u);

        // uがゼロベクトルの場合は変換をスキップ---------------
        if (w == 0.0) continue;
        
        // uの正規化-------------------------------------------
        for (j = i + 1; j < N; j++)  u[j] /= w;

        // p = A * u を計算------------------------------------
        matrix_vector_product(p, a, u);

        // u^t * p を計算--------------------------------------
        w = inner_product(u, p);
    
        // q = 2 * (p - inner * u) を計算----------------------
        for (j = 0; j < N; j++) {
            q[j] = 2.0 * (p[j] - w * u[j]);
        }

        // A = A - u * q^t - q * u^t を計算--------------------
        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                a[j][k] -= u[j] * q[k] + q[j] * u[k];
            }
        }
    }
}

// メイン関数-----------------------------------------------------------
int main(){

    // 対称行列A（鏡映変換の対象）
    double a[N][N] = {
        {  2.0,  6.0,  0.0, -2.0,  4.0 },
        {  6.0, -2.0, -1.0,  3.0, -7.0 },
        {  0.0, -1.0,  4.0,  9.0, -3.0 },
        { -2.0,  3.0,  9.0,  1.0, -5.0 },
        {  4.0, -7.0, -3.0, -5.0,  3.0 }
    };

    householder_trans(a);

    // 結果出力------------------------------------------------
    printf("三重対角化後の行列：\n");
    int j, k;
    for (j = 0; j < N; j++) {
        for (k = 0; k < N; k++) {
            printf("%12.6lf ", a[j][k]);
        }
        printf("\n");
    }

    return 0;
}
