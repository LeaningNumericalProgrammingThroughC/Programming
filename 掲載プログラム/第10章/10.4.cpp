#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 4           // ������
#define EPS 1e-8      // ��������p臒l

// ���όv�Z---------------------------------------------------------------------
double inner_product(double x[], double y[]) {

    double s = 0.0;
    int i;

    for (i = 0; i < N; i++) {
        s += x[i] * y[i];
    }

    return s;
}

// �x�N�g����2�m����------------------------------------------------------------
double vector_norm2(double x[]) {

    double s = 0.0;
    int i;
    
    for (i = 0; i < N; i++) {
        s += x[i] * x[i];
    }

    return sqrt(s);
}

// LU�����i�����s�{�b�g�I��t���j-----------------------------------------------
void lu_decomposition(double a[][N], int p[]) {

    double m, max, tmp;
    int i, j, k, max_i;

    for (i = 0; i < N - 1; i++) {

        max = fabs(a[i][i]);
        max_i = i;
        for (j = i + 1; j < N; j++) {
            if (fabs(a[j][i]) > max) {
                max = fabs(a[j][i]);
                max_i = j;
            }
        }
        p[i] = max_i;

        if (max < EPS) {
            printf("�W���s�񂪐����ł͂���܂���\n");
            exit(1);
        }
        
        // �s�̓��ꊷ��
        if (max_i != i) {
            for (j = 0; j < N; j++) {
                tmp = a[i][j]; a[i][j] = a[max_i][j]; a[max_i][j] = tmp;
            }
        }
        
        // �O�i����
        for (j = i + 1; j < N; j++) {
            m = a[j][i] / a[i][i];
            for (k = i + 1; k < N; k++) {
                a[j][k] -= m * a[i][k];
            }
            a[j][i] = m;
        }
    }
}

// �O�i����ƌ�ޑ���ɂ��A���������̉�@-------------------------------------
void lu_substitution(double a[][N], double b[], int p[]) {

    double tmp, s;
    int i, j;

    // �s�{�b�g�I���ɂ��b�̓���ւ�
    for (i = 0; i < N - 1; i++) {
        tmp = b[i]; b[i] = b[p[i]]; b[p[i]] = tmp;
    }

    // �O�i����i���O�p�s�� L �̉�@�j
    for (i = 0; i < N; i++) {

        for (j = 0; j < i; j++) {
            b[i] -= a[i][j] * b[j];
        }
    }

    // ��ޑ���i��O�p�s�� U �̉�@�j
    b[N - 1] = b[N - 1] / a[N - 1][N - 1];
    for (i = N - 2; i >= 0; i--) {

        s = 0.0;
        for (j = i + 1; j < N; j++) {
            s += a[i][j] * b[j];
        }

        b[i] = (b[i] - s) / a[i][i];
    }
}

// �t�����@------------------------------------------------------------------------
void inverse_iteration(const double A[][N], double lambda, double eigenvector[N],
    double* corrected_lambda) {

    double b[N][N], x[N], x_new[N], r = 0.0, r_old = 0.0,norm;
    int p[N], i, j;

    // b = A - lambda * I �̌v�Z
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            b[i][j] = A[i][j];
        }
        b[i][i] -= lambda;
    }

    // �ŗL�x�N�g���ɏ����l���Z�b�g
    x[0] = 1.0;
    for (i = 1; i < N; i++) {
        x[i] = 0.0;
    }

    // LU����
    lu_decomposition(b, p);

    do {
        
        for (i = 0; i < N; i++) x_new[i] = x[i];
        
        // �O�i����y�ь�ޑ��
        lu_substitution(b, x_new, p);

        r_old = r;
        r = inner_product(x_new, x);

        // x_new�𐳋K������
        norm = vector_norm2(x_new);
        for (i = 0; i < N; i++)  x[i] = x_new[i] / norm;
 
    } while (fabs((r - r_old) / r) > EPS);

    // �߂�l�̃Z�b�g-------------------------------------------------------
    *corrected_lambda = lambda + 1.0 / r;
    for (i = 0; i < N; i++) eigenvector[i] = x[i];
}

// ���C���֐�-----------------------------------------------------------------------
int main(void) {

    double A[N][N] = {
        {  3.0, -2.0,  1.0, -6.0 },
        {  3.0,  5.0, -4.0,  3.0 },
        {  7.0, -2.0, -3.0,  1.0 },
        { -4.0,  1.0,  3.0,  8.0 }
    };
    double lambda[N] = { 10.711, -5.823, 5.372, 2.740 };
    double eigenvector[N];
    double corrected_lambda;
    int i, k;

    printf("�ŗL�l : �ŗL�x�N�g��\n");
    for (k = 0; k < N; k++) {

        // �t�����@�����s---------------------------------------------------
        inverse_iteration(A, lambda[k], eigenvector, &corrected_lambda);
        
        // ���ʏo��---------------------------------------------------------
        printf("%10.6lf : [", corrected_lambda);
        for (i = 0; i < N - 1; i++) {
            printf("%10.6lf, ", eigenvector[i]);
        }
        printf("%10.6lf]\n", eigenvector[N - 1]);
    }

    return 0;
}
