#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 4          // ������ 
#define EPS 1e-10    // ����������p臒l 

// �����s�{�b�g�I��t���K�E�X�̏����@----------------------------------------------
void pivoted_gauss(double a[][N], double b[], double x[])
{
    int i, j, k, max_i;
    double m, s, max, tmp;

    // �O�i����------------------------------------------------------------
    for (i = 0; i < N - 1; i++)
    {
        // �����s�{�b�g�I��-----------------------------------------
        max = fabs(a[i][i]);
        max_i = i;
        for (j = i + 1; j < N; j++) {

            if (fabs(a[j][i]) > max) {
                max = fabs(a[j][i]);
                max_i = j;
            }
        }

        if (max < EPS) {

            printf("�W���s�񂪐����ł͂���܂���");
            exit(1);
        }
        else if (max_i != i) {

            for (j = i; j < N; j++) {
                tmp = a[i][j]; a[i][j] = a[max_i][j]; a[max_i][j] = tmp;
            }
            tmp = b[i]; b[i] = b[max_i]; b[max_i] = tmp;
        }

        // �O�i����-------------------------------------------------
        for (j = i + 1; j < N; j++) {

            m = a[j][i] / a[i][i];

            for (k = i + 1; k < N; k++) {
                a[j][k] = a[j][k] - m * a[i][k];
            }
            b[j] = b[j] - m * b[i];
        }
    }

    // ��ޑ��------------------------------------------------------------
    x[N - 1] = b[N - 1] / a[N - 1][N - 1];
    for (i = N - 2; i >= 0; i--) {

        s = 0;
        for (j = i + 1; j < N; j++) {
            s += a[i][j] * x[j];
        }
        x[i] = (b[i] - s) / a[i][i];
    }
}

// ���C���֐�----------------------------------------------------------------------
int main(void){

    double a[N][N] = {
        {  0.0, -3.0, -2.0, -4.0 },
        {  3.0, -2.0,  2.0, -1.0 },
        { -2.0,  0.0,  3.0,  2.0 },
        {  5.0,  2.0, -3.0,  3.0 }
    };
    double b[N] = { 2.0, -3.0, -2.0, 1.0 }, x[N];
    int i;

    pivoted_gauss(a, b, x);

    printf("���͎��̒ʂ�\n");
    for (i = 0; i < N; i++) {
        printf("x[%d] = %10.6lf\n", i, x[i]);
    }

    return 0;
}

