#include <stdio.h>
#include <math.h>

#define N 4  // �s��̎���

// �s���1�m����--------------------------------------------
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

// �s��̍ő�l�m����---------------------------------------
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

// ���C���֐�-----------------------------------------------
int main(){

    double a[N][N] = {
        {  1.0, -2.0,  3.0, -4.0},
        {  0.0,  1.0,  2.0, -1.0},
        { -1.0, -3.0,  4.0,  1.0},
        {  2.0,  0.0,  3.0,  5.0}
    };

    printf("1�m����      : %lf\n", matrix_norm1( a ));
    printf("������m���� : %lf\n", matrix_norm_max( a ));

    return 0;
}
