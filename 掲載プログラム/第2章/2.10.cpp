#include <stdio.h>
#include <math.h>

#define N 4  // s—ñ‚ÌŸ”

// s—ñ‚Ì1ƒmƒ‹ƒ€--------------------------------------------
double matrix_norm1(double a[][N]){

    double s, norm;
    int i, j;

    norm = 0.0;
    for (i = 0; i < N; i++) {

        s = 0.0;
        for (j = 0; j < N; j++) {
            s += fabs(a[j][i]);
        }
        
        if (s > norm) {
            norm = s;
        }
    }

    return norm;
}

// s—ñ‚ÌÅ‘å’lƒmƒ‹ƒ€---------------------------------------
double matrix_norm_max(double a[][N]){

    double s, norm;
    int i, j;

    norm = 0.0;
    for (i = 0; i < N; i++) {
        
        s = 0.0;
        for (j = 0; j < N; j++) {
            s += fabs(a[i][j]);
        }

        if (s > norm) {
            norm = s;
        }
    }

    return norm;
}

// ƒƒCƒ“ŠÖ”-----------------------------------------------
int main(){

    double a[N][N] = {
        {  1.0, -2.0,  3.0, -4.0},
        {  0.0,  1.0,  2.0, -1.0},
        { -1.0, -3.0,  4.0,  1.0},
        {  2.0,  0.0,  3.0,  5.0}
    };

    printf("1ƒmƒ‹ƒ€      : %lf\n", matrix_norm1( a ));
    printf("–³ŒÀ‘åƒmƒ‹ƒ€ : %lf\n", matrix_norm_max( a ));

    return 0;
}
