#include <stdio.h>
#include <math.h>

#define N 4

// 狭義優対角かチェックする関数---------------------------------------
int is_sdd(double a[N][N]) {
    
    double d, s;
    int i,j;

    for (i = 0; i < N; i++) {
        
        d = fabs(a[i][i]);
        
        s = 0.0;
        for (j = 0; j < N; j++) {
            if (i != j) s += fabs(a[i][j]);
        }

        if (d <= s) return 0;
    }

    return 1;
}

// メイン関数---------------------------------------------------------
int main() {
    double a[N][N] = {
        {  5.0, -1.0,  2.0, -1.0 },
        { -2.0,  7.0, -1.0,  3.0 },
        {  1.0,  2.0,  4.0,  3.0 },
        {  1.0,  4.0, -2.0,  9.0 },
    };

    if (is_sdd(a)) {
        printf("狭義優対角行列です．\n");
    }
    else {
        printf("狭義優対角行列ではありません．\n");
    }

    return 0;
}
