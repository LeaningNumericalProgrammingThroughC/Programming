#include <stdio.h>
#include <math.h>

#define N 5

// 1�m����-----------------------------------------------------
double vector_norm1( double x[] ){

    double s;
    int i;

    s = 0.0;
    for (i = 0; i < N; i++) {
        s += fabs(x[i]);
    }

    return s;
}

// 2�m����-----------------------------------------------------
double vector_norm2(double x[]){
    
    double s;
    int i;

    s = 0.0;
    for (i = 0; i < N; i++) {
        s += x[i] * x[i];
    }

    return sqrt(s);
}

// �ő�l�m����------------------------------------------------
double vector_norm_max(double x[]){

    double max;
    int i;

    max = fabs(x[0]);
    for (i = 1; i < N; i++) {
        if (fabs(x[i]) > max) {
            max = fabs(x[i]);
        }
    }
    
    return max;
}

// ���C���֐�--------------------------------------------------
int main(void){

    double x[N] = { 4.5, -3.1, 8.2, 3.3, 1.7 };

    printf("1�m�����̒l      : %lf\n", vector_norm1(x));
    printf("2�m�����̒l      : %lf\n", vector_norm2(x));
    printf("�ő�l�m�����̒l : %lf\n", vector_norm_max(x));

    return 0;
}
