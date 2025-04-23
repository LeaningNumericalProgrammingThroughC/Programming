#include <stdio.h>
#include <math.h>

#define N    4          // ������
#define EPS    1e-10    // ��������p
#define MAX_ITER  500   // �ő唽����

// �s��ƃx�N�g���̏�Z-----------------------------------------------
void mat_vec_product(double r[], double a[][N], double x[]) {
    int i, j;

    for (i = 0; i < N; i++) {
        r[i] = 0.0;
        for (j = 0; j < N; j++) {
            r[i] += a[i][j] * x[j];
        }
    }
}

// �x�N�g���̓���-----------------------------------------------------
double inner_product(double a[N], double b[N]) {
    int i;
    double s = 0.0;

    for (i = 0; i < N; i++) {
        s += a[i] * b[i];
    }

    return s;
}

// �ׂ���@-----------------------------------------------------------
double power_method(double A[][N], double x[], int* r) {

    double Ax[N], lambda, x_sq, norm_x;
    int i, iter = 0;

    do {

        // �s��A�ƃx�N�g��x�̐ς��v�Z����
        mat_vec_product(Ax, A, x);
        
        // ���C���[���ɂ��ߎ��ŗL�l�����߂�
        lambda = inner_product(Ax, x);
        
        // Ax��2�m�������v�Z���� 
        x_sq = inner_product(Ax, Ax);
        norm_x = sqrt(x_sq);

        // ���̔����̂��߂ɁCAx�̒l��x�ɖ߂�
        for (i = 0; i < N; i++) x[i] = Ax[i] / norm_x;

        iter++;

    } while (fabs(x_sq - lambda*lambda) >= EPS && iter < MAX_ITER);

    *r = iter;
    return lambda;
}

// ���C���֐�---------------------------------------------------------
int main() {

    double A[N][N] = {
        {  8.0,  3.0, -1.0, -3.0},
        {  2.0, -1.0,  3.0,  4.0},
        {  2.0, -6.0,  2.0,  1.0},
        { -1.0, -2.0,  7.0, -2.0}
    };
    double x[N] = { 1.0, 1.0, 1.0, 1.0 }, lambda;
    int i, iter;

    lambda = power_method(A, x, &iter);

    // ���ʂ̕\��------------------------------------------
    if (iter == MAX_ITER) {
        printf("�ׂ���@���������܂���ł����D");
    } else {
        printf("�����񐔂�%d��ł��D\n", iter);
        printf("�ő�ŗL�l�� %f �ł��D\n", lambda);
        printf("�ŗL�x�N�g���͈ȉ��ł��C\n");
        for (i = 0; i < N; i++) {
            printf("x[%d] = %f\n",i, x[i]);
        }
    }
}
