#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N   4          // 次元数
#define EPS   1e-8     // 収束判定用閾値

// 内積計算-------------------------------------------------------------------
double inner_product(const double x[], const double y[]) {

    double s = 0.0;
    for (int i = 0; i < N; i++) {
        s += x[i] * y[i];
    }

    return s;
}

// ベクトルの2ノルム----------------------------------------------------------
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

// 行列のコピー---------------------------------------------------------------
void matrix_copy(double r[][N], double x[][N]) {

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            r[i][j] = x[i][j];
        }
    }
}

// 符号関数（ただし，0の場合は1を返す）---------------------------------------
int sign(double x) {
    return (x >= 0.0) ? 1 : -1;
}

// ハウスホルダー変換による三重対角化-----------------------------------------
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

// N(c)の計算-------------------------------------------------------------
int n(const double alpha[], double beta[], double c) {

    double w, p_k, p_k_2 = 1.0, p_k_1 = c - alpha[0];
    int i, count = 0;

    if (p_k_1 < 0) count++;
    if (p_k_1 != 0.0) w = p_k_1;

    for (i = 1; i < N; i++) {

        p_k = (c - alpha[i]) * p_k_1 - beta[i - 1] * beta[i - 1] * p_k_2;

        if (p_k * w < 0) count++;
        if (p_k != 0.0) w = p_k;

        p_k_2 = p_k_1;
        p_k_1 = p_k;
    }

    return count;
}

// 2分法------------------------------------------------------------------
double bisection_eigenvalue(int index, double alpha[], double beta[],
    double* min, double* max) {

    double c, lower = *min, upper = *max;

    while (upper - lower > EPS) {

        c = (upper + lower) / 2.0;
        if (n(alpha, beta, c) >= index) {
            lower = c;
        }
        else {
            upper = c;
        }
    }

    /* 戻り値をセット */
    *min = lower; *max = upper;
    return (upper + lower) / 2.0;
}

// 2分法による全固有値の計算--------------------------------------------------
void eigenvalues(double lambda[], double a[][N]) {

    double alpha[N], beta[N - 1];
    double min_0, min, max;
    int i;

    for (i = 0; i < N; i++) alpha[i] = a[i][i];
    for (i = 0; i < N - 1; i++) beta[i] = a[i][i + 1];
    double c, w;

    // 初期区間の設定----------------------------------------------
    min = alpha[0] - fabs(beta[0]);
    max = alpha[0] + fabs(beta[0]);

    for (i = 1; i < N; i++) {

        if (i == N - 1) {
            w = fabs(beta[i - 1]);
        } else {
            w = fabs(beta[i - 1]) + fabs(beta[i]);
        }

        if ((alpha[i] - w) < min) {
            min = alpha[i] - w;
        }
        if ((alpha[i] + w) > max) {
            max = alpha[i] + w;
        }
    }
    min_0 = min;

    // 全ての固有値を2分法で求める----------------------------------
    for (i = 1; i <= N; i++) {

        lambda[i - 1] = bisection_eigenvalue(i, alpha, beta,
            &min, &max);
        max = min; min = min_0;
    }
}

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
    for (i = 0; i < N - 1; i++) {
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

// 逆反復法-------------------------------------------------------------------
void inverse_iteration_method(double a[][N], double lambda[]) {

    double x_k[N] = { 1,1,1,1 }, x_k_1[N], b[N][N];
    double w, r = 0.0, r_old;
    int p[N], i, j, k;

    printf("固有値：固有ベクトル\n");

    for (k = 0; k < N; k++) {

        // A - λＩの計算---------------------------------------
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) b[i][j] = a[i][j];
            b[i][i] -= lambda[k];
        }

        // 固有ベクトルの初期値設定-----------------------------
        x_k[0] = 1.0;
        for (i = 1; i < N; i++) x_k[i] = 0.0;

        // LU分解
        lu_decomposition(b, p);

        do {

            // LU分解によりx_k+1を求める-------------------------
            for (i = 0; i < N; i++) x_k_1[i] = x_k[i];
            lu_substitution(b, x_k_1, p);

            // レイリー商の計算----------------------------------
            r_old = r;
            r = inner_product(x_k_1, x_k);

            // 固有ベクトルの正規化------------------------------
            w = vector_norm2(x_k_1);
            for (i = 0; i < N; i++) x_k[i] = x_k_1[i] / w;

        } while (fabs((r - r_old) / r) > EPS);

        // 結果出力-----------------------------------------------
        printf("%12.8lf : [", 1. / r + lambda[k]);
        for (i = 0; i < N - 1; i++) {
            printf("%12.8lf,", x_k[i]);
        }
        printf("%12.8lf ]^t\n", x_k[N - 1]);
    }
}

// メイン関数-----------------------------------------------------------
int main() {

    // 対称行列A（鏡映変換の対象）
    double a[N][N] = {
        {  3.0,  2.0,  0.0,  1.0 },
        {  2.0,  5.0,  2.0,  7.0 },
        {  0.0,  2.0,  2.0,  6.0 },
        {  1.0,  7.0,  6.0,  1.0 },
    };
    double a_org[N][N], lambda[N];
    int i, j;

    matrix_copy(a_org, a);

    // ハウスホルダー法
    householder_trans(a);

    // 3重対角行列を出力
    printf("3重対角化後の行列：\n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%12.6lf ", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    // 2分法による固有値計算
    eigenvalues(lambda, a);

    // 逆反復法
    inverse_iteration_method(a_org, lambda);

    return 0;
}

