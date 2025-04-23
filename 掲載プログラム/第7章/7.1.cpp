#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define EPS  1e-8         /* 収束判定用 */
#define MAX_ITER  100     /* 反復回数上限 */

// 関数の定義----------------------------------------------------------
double f(double x) {
    return  pow(x, 6.0) + 3 * pow(x, 4.0) - pow(x, 3.0) - 7.0;
}

// 2分法---------------------------------------------------------------
double bisection(double a, double b, int* r) {

    double c, f_a, f_b, f_c;
    int iter = 0;
    
    // 同じ関数値の再計算を避けるため変数に保存して活用する
    f_a = f(a); f_b = f(b); 

    // 初期チェック-------------------------------------------
    if (f_a == 0.0) { *r = iter; return a; }
    if (f_b == 0.0) { *r = iter; return b; }
    if (f_a * f_b > 0.0) {
        printf("解の存在を判定できません．\n");
        exit(0); 
    }

    // 反復改良-----------------------------------------------
    while (fabs(b - a) > EPS) {
        
        iter++;
        c = (a + b) / 2.0;
        f_c = f(c);

        if (f_c == 0.0) break;

        if (f_a * f_c < 0.0) {
            b = c;
        } else {
            a = c; 
            f_a = f_c;
        }
    }

    *r = iter;
    return (a + b) / 2.0;
}

// メイン関数----------------------------------------------------------
int main() {

    double a = 0.0, b = 5.0; // 初期区間[a,b]の設定
    int iter;

    double c = bisection(a, b, &iter);

    printf("反復回数は %d 回，解は %lf です．\n", iter, c);
    
    return 0;
}

