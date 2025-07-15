#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <cmath>

#define N      5            // ルジャンドル多項式の次数
#define NP   (N+1)          // 標本点数 = N+1
#define EPS_BISECTION 1e-1  // 2分法用の収束判定閾値
#define EPS_DKA 1e-8        // DKA法用の収束判定閾値
#define MAX_ITER 100        // 最大反復回数

// 階乗計算（戻り値:double）-------------------------------------------------------
double factorial(int n) {

    if (n < 0) return 0.0;

    unsigned long long p, i;

    p = 1;
    for (i = 2; i <= n; i++) p *= i;

    return (double)p;
}

// 複素数クラス--------------------------------------------------------------------
class Complex {

    // データメンバ------------------------------------------------------
    double Re, Im;

public:
    /* コンストラクタ */
    Complex(double x = 0, double y = 0) : Re(x), Im(y) { }

    /* 出力関数 */
    void print() {
        printf("%.6lf + %.6lfi ", Re, Im);
    }

    /* データへのアクセス */
    double get_Re() const { return Re; }
    double get_Im() const { return Im; }

    /* 絶対値の計算 */
    double abs() const {
        return sqrt(Re * Re + Im * Im);
    }

    // 演算子オーバーロード---------------------------------------------
    /* マイナス演算子 */
    Complex operator-() const {
        return Complex(-Re, -Im);
    }
    /* 加算 */
    friend Complex operator+(const Complex& x, const Complex& y) {
        return Complex(x.Re + y.Re, x.Im + y.Im);
    }
    /* 減算 */
    friend Complex operator-(const Complex& x, const Complex& y) {
        return Complex(x.Re - y.Re, x.Im - y.Im);
    }
    /* 乗算 */
    friend Complex operator*(const Complex& x, const Complex& y) {
        return Complex(x.Re * y.Re - x.Im * y.Im,
            x.Re * y.Im + x.Im * y.Re);
    }
    /* 除算 */
    friend Complex operator/(const Complex& x, const Complex& y) {
        double work = y.Re * y.Re + y.Im * y.Im;

        return Complex((x.Re * y.Re + x.Im * y.Im) / work,
            (-x.Re * y.Im + x.Im * y.Re) / work);
    }
    /* 指数関数 */
    friend Complex exp(const Complex& x) {
        double eR = std::exp(x.Re);
        return Complex(eR * cos(x.Im), eR * sin(x.Im));
    }
};

// 多項式P(z)の計算----------------------------------------------------------------
Complex P(double a[], Complex x) {
    int i;
    Complex r;

    r = a[0];
    for (i = 1; i <= NP; i++) {
        r = r * x + a[i];
    }
    return r;
}

// DK法の修正量delta_xの分母の計算-------------------------------------------------
Complex PI(Complex x[], int i) {
    int j;
    Complex r(1, 0);

    for (j = 0; j < NP; j++) {
        if (i != j) {
            r = r * (x[i] - x[j]);
        }
    }
    return r;
}

// 多項式Sの計算-------------------------------------------------------------------
double S(double b[], double x) {
    int i;
    double r;

    r = b[0];
    for (i = 1; i <= NP; i++) {
        r = r * x + b[i];
    }
    return r;
}

// 2分法---------------------------------------------------------------------------
double bisection_DKA(double b[], double left, double right) {

    double center, S_center;
    double S_left = S(b, left);

    while ((right - left) > EPS_BISECTION * right) {

        center = (right + left) / 2.0;
        S_center = S(b, center);

        if (S_center == 0.0) break;

        if (S_left * S_center < 0.0) {
            right = center;
        }
        else {
            left = center;
            S_left = S_center;
        }
    }

    return (left + right) / 2.0;
}

// DKA法---------------------------------------------------------------------------
int dka(double a[], Complex* z, int* r) {

    Complex delta_z;
    double R, b[NP + 1], max_delta_z, zeta, eta = 0, tmp;
    int i, j, iter = 0, converged = 0;

    /* 解の重心を計算 */
    zeta = -a[1] / (a[0] * NP);

    for (i = 0; i < NP + 1; i++)  b[i] = a[i];

    /* 多項式P^*の係数を計算 */
    for (i = 0; i <= NP - 2; i++) {
        for (j = 1; j <= NP - i; j++) {
            b[j] = b[j] + b[j - 1] * zeta;
        }
    }

    /* 多項式Sの係数を計算 */
    int m = 0;
    for (i = 2; i <= NP; i++) {
        b[i] = -fabs(b[i]);
        if (fabs(b[i]) > 1.0e-8) m++;
    }

    /* ηの計算 */
    for (i = 2; i <= NP; i++) {
        tmp = pow(m * fabs(b[i]) / fabs(b[0]), 1.0 / (double)i);
        if (tmp > eta) eta = tmp;
    }

    /* 2分法の実行 */
    R = bisection_DKA(b, 0, eta);

    /* 初期値の設定 */
    for (j = 0; j < NP; j++) {
        z[j] = zeta + R * exp(Complex(0, 1) * (2 * M_PI * ((double)j) / NP)
            + M_PI / (2 * NP));
    }

    /* DK法の反復更新 */
    while (iter < MAX_ITER) {
        max_delta_z = 0;
        for (i = 0; i < NP; i++) {

            delta_z = P(a, z[i]) / (a[0] * PI(z, i));
            z[i] = z[i] - delta_z;

            if (delta_z.abs() > max_delta_z) {
                max_delta_z = delta_z.abs();
            }
        }
        iter++;

        // 収束判定---------------------------
        if (max_delta_z <= EPS_DKA) {
            converged = 1;
            break;
        }
    }

    /* 戻り値をセット */
    *r = iter;
    return converged;
}

// 重みwの計算-------------------------------------------------------------------------------
double w(int n, double xi, double Pn1_at_xi) {

    return 2.0 * (1.0 - xi * xi) / ((double)n * (double)n * Pn1_at_xi * Pn1_at_xi);
}

// ルジャンドル多項式の値を計算--------------------------------------------------------------
double P_n_x(int n, double x) {

    double P_k_m_1, P_k, P_k_p_1;

    if (n == 0) return 1.0;
    if (n == 1) return x;

    P_k_m_1 = 1.0;
    P_k = x;

    for (int j = 1; j <= n - 1; j++) {

        P_k_p_1 = (2.0 * (double)j + 1.0) / ((double)j + 1.0) * x * P_k
            - (double)j / ((double)j + 1.0) * P_k_m_1;

        P_k_m_1 = P_k;
        P_k = P_k_p_1;
    }

    return P_k_p_1;
}

// 関数の定義-----------------------------------------------------------
double func(double x) {
    
    return x + cos(x);
}

// メイン関数-----------------------------------------------------------
int main() {

    Complex z[NP];
    double x_nodes[NP], s, c[NP + 1];
    int t, i, iter, result;

    // N+1次のルジャンドル多項式の係数c[i]の値をセットする--------------------------------
    for (i = 0; i <= NP; i++) {

        if ((NP - i) % 2 == 1) {    
            c[NP - i] = 0.0;
        } else {
            t = (NP - i) / 2;
            c[NP - i] = pow(-1, t) * factorial(2 * (NP - t))
                / (pow(2, NP) * factorial(t) * factorial(NP - t) * factorial(i));
        }
    }

    // N+1次のルジャンドル多項式の零点をDKA法で求める-------------------------------------
    result = dka(c, z, &iter);

    // DKA法が収束しなかった場合は処理を終了する------------------------------------------
    if (result == 0) {
        printf("DKA法の反復回数が上限に達しました．\n");
        exit(0);
    }

    // DKA法の解の実部を零点とする（虚部は0）---------------------------------------------
    for (i = 0; i < NP; i++) {
        x_nodes[i] = z[i].get_Re();
    }

    // 積分値の計算（式(13.22））---------------------------------------------------------
    s = 0.0;
    for (i = 0; i < NP; i++) {

        s += w(NP, x_nodes[i], P_n_x(N, x_nodes[i])) * func(x_nodes[i]);
    }

    // ガウス・ルジャンドル公式
    printf("ガウス・ルジャンドル公式 ： %12.10lf\n", s);

    return 0;
}
