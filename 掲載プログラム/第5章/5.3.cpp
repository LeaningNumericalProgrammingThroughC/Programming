#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N   4           /* ������ */
#define EPS  1e-8       /* epsilon�̐ݒ� */
#define MAX_ITER 10     /* �ő唽���� */

// ���όv�Z-------------------------------------------------------------------
double inner_product(const double x[], const double y[]) {

    double s = 0.0;
    for (int i = 0; i < N; i++) {
        s += x[i] * y[i];
    }

    return s;
}

// �ő�l�m����---------------------------------------------------------------
double vector_norm_max(double x[]) {

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

// �s��ƃx�N�g���̐�---------------------------------------------------------
void matrix_vector_product(double c[], double a[][N], double b[]) {

    double s;
    int i, j;

    for (i = 0; i < N; i++){

        s = 0.0;
        for (j = 0; j < N; j++){

            s += a[i][j] * b[j];
        }
        c[i] = s;
    }
}

// �������z�@-----------------------------------------------------------------
void cg(double a[][N], double b[], double x[]) {

    double r[N], p[N], alpha, beta, tmp1[N], tmp2;
    int i, iter = 1;

    // x^(1) �̌v�Z------------------------------------------------------
    matrix_vector_product(tmp1, a, x);
    for (i = 0; i < N; i++) {
        r[i] = b[i] - tmp1[i];
        p[i] = r[i];
    }

    // A * p_k �̌v�Z 
    matrix_vector_product(tmp1, a, p);
    // (p_k, A*p_k) �̌v�Z 
    tmp2 = inner_product(p, tmp1);

    alpha = inner_product(p, p) / tmp2;

    for (i = 0; i < N; i++) x[i] = x[i] + alpha * p[i];
    for (i = 0; i < N; i++) r[i] = r[i] - alpha * tmp1[i];

    // x^(2) �ȍ~�̔����v�Z----------------------------------------------
    while (iter++ < MAX_ITER){

        beta = -inner_product(r, tmp1) / tmp2;
        for (i = 0; i < N; i++) p[i] = r[i] + beta * p[i];

        // A * p_k �̌v�Z 
        matrix_vector_product(tmp1, a, p);
        // (p_k, A*p_k) �̌v�Z 
        tmp2 = inner_product(p, tmp1);

        alpha = inner_product(p, r) / tmp2;

        for (i = 0; i < N; i++) x[i] = x[i] + alpha * p[i];
        for (i = 0; i < N; i++) r[i] = r[i] - alpha * tmp1[i];

        // ��������--------------------------------------------------
        if (vector_norm_max(r) <= EPS) {
            break;
        }
    }

    // ���ʕ\��----------------------------------------------------------
    if (iter > MAX_ITER) {
        printf("�������܂���ł���\n");
    }
    else {
        printf("���͈ȉ��̒ʂ�\n");
        printf("�����񐔁F%d��\n", iter);
        for (i = 0; i < N; i++) {
            printf("x_k[%d] = %lf\n", i, x[i]);
        }
    }
}

// ���C���֐�-----------------------------------------------------------------
int main(void){

    double a[N][N] = {
        { 6.0, 2.0, 1.0, 1.0 },
        { 2.0, 7.0, 3.0, 3.0 },
        { 1.0, 3.0, 8.0, 2.0 },
        { 1.0, 3.0, 2.0, 7.0 },
    };
    double b[N] = { 7.0, 6.0, -5.0, -3.0 };
    double x[N] = { 1.0, 1.0, 1.0, 1.0 };

    cg( a, b, x );

    return 0;
}



