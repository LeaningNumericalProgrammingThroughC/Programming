#include <stdio.h>
#include <stdint.h>
#include <math.h>

// メモリ上のビット列を出力
void print_bits(double x) {
    union { double d; uint64_t u; } a;
    a.d = x;

    // 見出しを出力
    printf(" s |       e       |                              f\n");
    printf("---+---------------+------------------------------------------------------------------\n");

    // 1ビットずる画面出力
    putchar(' ');
    for (int i = 63; i >= 0; i--) {

        // ビットの値を出力
        putchar((a.u & (1ULL << i)) ? '1' : '0');
        
        // 区切り文字を出力
        if (i == 63 || i == 52) printf(" | ");
        else if (i % 4 == 0 && i != 0) putchar(' ');
    }
    putchar('\n');
}

// メイン関数
int main(void) {

    union {
        uint64_t u;
        double   d;
    } x;
    
    // 無限大のビット列を出力
    double y = 0.0;
    x.d = 1.0 / y;

    printf("x = %.20e\n", x.d);
    print_bits(x.d);
    printf("\n");

    // 非正規化数の適当な値を一例としてビット列を出力
    x.d = pow(2, -1050);

    printf("x = %.20e\n", x.d);
    print_bits(x.d);

    return 0;
}
