#include <stdio.h>

#define N 11    // �W�{�_��

// ���O�����W�����---------------------------------------------------------
double lagrange(double sample_x[], double sample_y[],int n, double x) {

    double s = 0.0, p;
    int i, j;

    // �� f(x)*l(x)�̃��[�v-------------------------------------------
    for (i = 0; i < n; i++) {
        
        // ��(x-x_j)/(x_i-x_j)�̃��[�v--------------------------------
        p = sample_y[i];
        for (j = 0; j < i; j++) {
            p *= (x - sample_x[j]) / (sample_x[i] - sample_x[j]);
        }
        for (j = i+1; j < n; j++) {
            p *= (x - sample_x[j]) / (sample_x[i] - sample_x[j]);
        }
        s += p;
    }

    return s;
}

// ���C���֐�---------------------------------------------------------------
int main(void) {

    // �W�{�_�̒�`
    double sample_x[N] = { 1.0, 2.7, 3.0, 4.5, 5.2,
                                6.1, 7.4, 8.2, 9.3, 10.0, 11.5 };
    double sample_y[N] = { 3.5, 5.9, 6.5, 5.0, 4.1, 
                                3.2, 3.0, 4.3, 9.5, 8.3, 8.2 };
    double x, y, m = 100, k;
    double h = (sample_x[N - 1] - sample_x[0]) / m;

    // ��Ԃƌ��ʏo��------------------------------------------------
    printf("    x               y\n");
    printf("----------------------------------\n");
    for (k = 0; k <= m; k++) {
        x = sample_x[0] + k * h;
        y = lagrange(sample_x, sample_y, N, x);
        printf("%10.6lf\t%10.6lf\n", x, y);
    }

    return 0;
}
