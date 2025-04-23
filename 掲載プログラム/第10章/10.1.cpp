#include <stdio.h>
#include <math.h>

#define   N   5     // ������

// ���όv�Z------------------------------------------------------------------
double inner_product(const double x[], const double y[]) {

    double s = 0.0;
    for (int i = 0; i < N; i++) {
        s += x[i] * y[i];
    }

    return s;
}

// 2�m����-------------------------------------------------------------------
double vector_norm2(double x[]) {

    double s;
    int i;

    s = 0.0;
    for (i = 0; i < N; i++) {
        s += x[i] * x[i];
    }

    return sqrt(s);
}

// �s��ƃx�N�g���̐�---------------------------------------------------------
void matrix_vector_product(double c[], double a[][N], double b[]) {

    double s;
    int i, j;

    for (i = 0; i < N; i++) {

        s = 0.0;
        for (j = 0; j < N; j++) {

            s += a[i][j] * b[j];
        }
        c[i] = s;
    }
}

// �����֐��i�������C0�̏ꍇ��1��Ԃ��j---------------------------------
int sign(double x) {
    return (x >= 0.0) ? 1 : -1;
}

// ���f�ϊ��ɂ��O�d�Ίp��---------------------------------------------
void householder_trans(double a[][N]) {

    double s, w, alpha, u[N], p[N], q[N];
    int i, j, k;

    // i: ���˕ϊ���K�p�����i0�`N-3�܂Łj
    for (i = 0; i < N - 2; i++) {

        // u�̍\��---------------------------------------------
        s = 0.0;
        for (j = i + 1; j < N; j++) {
            s += a[j][i] * a[j][i];
        }
        alpha = -sign(a[i + 1][i]) * sqrt(s);

        for (j = 0; j <= i; j++) u[j] = 0.0;
        u[i + 1] = a[i + 1][i] - alpha;
        for (j = i + 2; j < N; j++) {
            u[j] = a[j][i];
        }

        w = vector_norm2(u);

        // u���[���x�N�g���̏ꍇ�͕ϊ����X�L�b�v---------------
        if (w == 0.0) continue;
        
        // u�̐��K��-------------------------------------------
        for (j = i + 1; j < N; j++)  u[j] /= w;

        // p = A * u ���v�Z------------------------------------
        matrix_vector_product(p, a, u);

        // u^t * p ���v�Z--------------------------------------
        w = inner_product(u, p);
    
        // q = 2 * (p - inner * u) ���v�Z----------------------
        for (j = 0; j < N; j++) {
            q[j] = 2.0 * (p[j] - w * u[j]);
        }

        // A = A - u * q^t - q * u^t ���v�Z--------------------
        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                a[j][k] -= u[j] * q[k] + q[j] * u[k];
            }
        }
    }
}

// ���C���֐�-----------------------------------------------------------
int main(){

    // �Ώ̍s��A�i���f�ϊ��̑Ώہj
    double a[N][N] = {
        {  2.0,  6.0,  0.0, -2.0,  4.0 },
        {  6.0, -2.0, -1.0,  3.0, -7.0 },
        {  0.0, -1.0,  4.0,  9.0, -3.0 },
        { -2.0,  3.0,  9.0,  1.0, -5.0 },
        {  4.0, -7.0, -3.0, -5.0,  3.0 }
    };

    householder_trans(a);

    // ���ʏo��------------------------------------------------
    printf("�O�d�Ίp����̍s��F\n");
    int j, k;
    for (j = 0; j < N; j++) {
        for (k = 0; k < N; k++) {
            printf("%12.6lf ", a[j][k]);
        }
        printf("\n");
    }

    return 0;
}
