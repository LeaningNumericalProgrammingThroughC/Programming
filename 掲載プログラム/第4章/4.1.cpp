#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N    4       // ������ 
#define EPS 1e-10    // ����������p

// LU����----------------------------------------------------------------------
void lu_decomposition(double a[][N], int p[]) {

    double m, s, max, tmp;
    int i, j, k, max_i;

    for (i = 0; i < N - 1; i++)
    {
        // �����s�{�b�g�I��-----------------------------------
        // �ő�s�̌���
        max = fabs(a[i][i]);
        max_i = i;
        for (j = i + 1; j < N; j++) {

            if (fabs(a[j][i]) > max) {
                max = fabs(a[j][i]);
                max_i = j;
            }
        }

        // ���ꊷ�����s��z��ɋL��
        p[i] = max_i;
        
        // �������̔���
        if (max < EPS) {
            printf("�W���s�񂪐����ł͂���܂���");
            exit(1);
        }
        
        // �s�̓��ꊷ��
        if (max_i != i) {

            for (j = 0; j < N; j++) {
                tmp = a[i][j]; a[i][j] = a[max_i][j]; a[max_i][j] = tmp;
            }
        }
      
        // �O�i����-------------------------------------------
        for (j = i + 1; j < N; j++) {

            m = a[j][i] / a[i][i];

            // �s��U�̗v�f���v�Z
            for (k = i + 1; k < N; k++) {
                a[j][k] = a[j][k] - m * a[i][k];
            }

            // �s��L�̗v�f���i�[
            a[j][i] = m;
        }
    }
}

// �O�i����ƌ�ޑ��----------------------------------------------------------
void lu_substitution(double a[][N], double b[], int p[]) {

    double tmp, s;
    int i, j, k;

    // �������̉E��Pb�����߂�-----------------------------------
    for (i = 0; i < N-1; i++) {
        tmp = b[i]; b[i] = b[p[i]]; b[p[i]] = tmp;
    }

    // �O�i���-------------------------------------------------
    for (i = 0; i < N; i++) {

        for (j = 0; j < i; j++) {
            b[i] -= a[i][j] * b[j];
        }
    }

    // ��ޑ��-------------------------------------------------
    b[N - 1] = b[N - 1] / a[N - 1][N - 1];
    for (i = N - 2; i >= 0; i--) {

        s = 0;
        for (j = i + 1; j < N; j++) {
            s += a[i][j] * b[j];
        }
        b[i] = (b[i] - s) / a[i][i];
    }
}

// ���C���֐�------------------------------------------------------------------
int main(void){

    double a[N][N] = {
        {  2.0,  2.0,  1.0,  2.0 },
        {  2.0,  3.0, -1.0,  3.0 },
        { -3.0,  4.0,  0.0,  3.0 },
        {  1.0,  3.0, -2.0,  1.0 }
    };
    double b[N] = { 7.0, 6.0, -5.0, -3.0 };
    int i, p[N];

    lu_decomposition(a, p);

    lu_substitution(a, b, p);

    // ���ʏo��-------------------------------------------------
    printf("���͎��̒ʂ�\n");
    for (i = 0; i < N; i++) {
        printf("x[%d] = %10.6lf\n", i, b[i]);
    }

    return 0;
}