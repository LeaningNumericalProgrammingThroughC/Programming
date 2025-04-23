#include <stdio.h>

#define N 4   /* ������ */

// �C���R���X�L�[����---------------------------------------------------------
void cholesky_decomp(double a[][N]) {

    double s;
    int i, j, k;

    for (i = 0; i < N; i++) {

        // �Ίp�s��D�����߂� (�s��A�̑Ίp�����Ɋi�[�j--------------------
        s = 0.0;
        for (k = 0; k < i; k++) {
            s += a[i][k] * a[k][k] * a[i][k];
        }
        a[i][i] -= s;

        // ���O�p�s��L�����߂�i�s��A�̉��O�p�����Ɋi�[�j----------------
        for (j = i + 1; j < N; j++){

            s = 0.0;
            for (k = 0; k < i; k++){

                s += a[i][k] * a[k][k] * a[j][k];
            }
            a[j][i] = (a[j][i] - s) / a[i][i];
        }
    }
}

// �O�i����ƌ�ޑ��---------------------------------------------------------
void cholesky_substitution(double a[][N], double b[N]) {

    double s;
    int i, j;

    // LD y = b�̉�y�����߂�iy�̒l��b�Ɋi�[�j------------------------
    for (i = 0; i < N; i++){

        s = 0.0;
        for (j = 0; j < i; j++){

            s += a[j][j] * a[i][j] * b[j];
        }
        b[i] = (b[i] - s) / a[i][i];
    }

    // L^t x = y �̉�x�����߂�ix�̒l��b�Ɋi�[)------------------------
    for (i = N - 2; i >= 0; i--) {

        s = 0.0;
        for (j = i + 1; j < N; j++){

            s += a[j][i] * b[j];
        }
        b[i] -= s;
    }
}

// ���C���֐�-----------------------------------------------------------------
int main(void) {

    double a[N][N] = {
        {  3.0, -2.0,  1.0,  2.0},
        { -2.0,  6.0,  0.0, -2.0},
        {  1.0,  0.0,  5.0,  4.0},
        {  2.0, -2.0,  4.0,  6.0}
    };
    double b[N] = { 7.0, -6.0, 9.0, 6.0 };
    int i;

    cholesky_decomp(a);

    cholesky_substitution(a, b);

    printf("���͎��̒ʂ�\n");
    for (i = 0; i < N; i++) {
        printf("x[%d] = %lf\n", i, b[i]);
    }

    return 0;
}