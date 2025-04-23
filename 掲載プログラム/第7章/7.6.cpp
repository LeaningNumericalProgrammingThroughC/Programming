#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES  // Visual Studio等でコンパイルする場合必要
#include <cmath>

#define N 5                    // 次数
#define EPS_BISECTION 1e-1     // 2分法用の収束判定閾値
#define EPS_DKA 1e-8           // DKA法用の収束判定閾値
#define MAX_ITER 100           // 最大反復回数


// 複素数クラス----------------------------------------------------------------
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

// 多項式P(z)の計算----------------------------------------------------------
Complex P(double a[], Complex x) {
    int i;
    Complex r;

    r = a[0];
    for (i = 1; i <= N; i++) {
        r = r * x + a[i];
    }
    return r;
}

// DK法の修正量delta_xの分母の計算-------------------------------------------
Complex PI(Complex x[], int i) {
    int j;
    Complex r(1, 0);

    for (j = 0; j < N; j++) {
        if (i != j) {
            r = r * (x[i] - x[j]);
        }
    }
    return r;
}

// 多項式Sの計算-------------------------------------------------------------
double S(double b[], double x) {
    int i;
    double r;

    r = b[0];
    for (i = 1; i <= N; i++) {
        r = r * x + b[i];
    }
    return r;
}

// 2分法---------------------------------------------------------------------
double bisection_DKA(double b[], double left, double right) {
    
    double center, S_center;
    double S_left = S(b, left);

    while ( (right - left) > EPS_BISECTION * right) {
        
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

// DKA法---------------------------------------------------------------------
int dka(double a[], Complex* z, int* r) {

    Complex delta_z;
    double R, b[N + 1], max_delta_z, zeta, eta = 0, tmp;
    int i, j, iter= 0, converged = 0;

    /* 解の重心を計算 */
    zeta = -a[1] / (a[0] * N);

    for (i = 0; i < N + 1; i++)  b[i] = a[i];

     /* 多項式P^*の係数を計算 */
    for (i = 0; i <= N-2; i++) {
        for (j = 1; j <= N - i; j++) {
            b[j] = b[j] + b[j - 1] * zeta;
        }
    }

    /* 多項式Sの係数を計算 */
    int m = 0;
    for (i = 2; i <= N; i++) {
        b[i] = -fabs(b[i]);
        if (fabs(b[i]) > 1.0e-8) m++;
    }

    /* ηの計算 */
    for (i = 2; i <= N; i++) {
        tmp = pow(m * fabs(b[i]) / fabs(b[0]), 1.0 / (double)i);
        if (tmp > eta) eta = tmp;
    }

    /* 2分法の実行 */
    R = bisection_DKA(b, 0, eta);
 
    /* 初期値の設定 */
    for (j = 0; j < N; j++) {
        z[j] = zeta + R * exp(Complex(0, 1) * (2 * M_PI * ((double) j) / N)
            + M_PI / (2 * N));
    }

    /* DK法の反復更新 */
    while (iter < MAX_ITER) {
        max_delta_z = 0;
        for (i = 0; i < N; i++) {

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

// メイン関数---------------------------------------------------------------
int main(void) {

    Complex z[N];
    // 多項式P(z)の係数（a[i] ： x^(N-i)の係数）
    double a[N + 1] = { 1,-2,4,-5,3,1 };
    int i,iter, result;

    result = dka(a, z, &iter);

    // 結果出力-----------------------------------------------------
    if (result == 1) {
        printf("解は以下の通りです．（反復回数：%d回）\n", iter);
        for (i = 0; i < N; i++) {
            printf("x[%d] = ", i); z[i].print(); printf("\n");
        }
    } else {
        printf("反復回数が上限に達しました．\n");
    }

    return 0;
}


