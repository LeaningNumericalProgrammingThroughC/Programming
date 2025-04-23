#include <stdio.h>

#define N 4  // ベクトルの次元数

// 内積計算-------------------------------------------------------
double inner_product(const double x[], const double y[]) {

    double s = 0.0;
    for (int i = 0; i < N; i++) {
        s += x[i] * y[i];
    }
    
    return s;
}

// メイン関数-----------------------------------------------------
int main(void) {

    double x[N] = { 1.0, 2.0, 3.0, 4.0 };
    double y[N] = { 5.0, 6.0, 7.0, 8.0 };

    double r = inner_product(x, y);

    printf("(x,y) = %lf\n", r);

    return 0;
}
