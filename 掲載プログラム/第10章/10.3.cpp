#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N    4          // �s��̎���
#define EPS   1e-8      // ��������p臒l
#define MAX_ITER 1000   // �ő唽����

// QR�@-----------------------------------------------------------------------------
void qr_method(double a[N][N]) {

    double Qt[N][N], sin_theta, cos_theta, tmp, w, work[N];
    int i, j, k, m = N - 1, iter = 0;

    while (1) {

        // �I������P�F�Ίp������1���̗v�f���S��0�ɂȂ�ΏI��-----------------
        while (m > 0 && fabs(a[m][m - 1]) <= EPS) m--;
        if (m == 0) break;

        // �I������Q�F�����񐔂ɂ��I������-----------------------------------
        if (iter >= MAX_ITER) {
            printf("�����񐔂�����ɒB���܂����D\n");
            exit(1);
        }

        // Q^t��P�ʍs��ŏ�����------------------------------------------------
        for (i = 0; i <= m; i++) {
            for (j = 0; j <= m; j++) Qt[i][j] = 0.0;
            Qt[i][i] = 1.0;
        }

        // QR�������s��---------------------------------------------------------
        for (i = 0; i <= m - 1; i++) {

            w = sqrt(a[i][i] * a[i][i] + a[i+1][i] * a[i+1][i]);
            
            /* w��0�̏ꍇ�́C�P�ʍs��Ƃ��� */
            if (fabs(w) < EPS) {
                sin_theta = 0.0;
                cos_theta = 1.0;
            } else {
                sin_theta = -a[i+1][i] / w;
                cos_theta = a[i][i] / w;
            }

            // R = P_i * R------------------------------------------------
            a[i][i] = w; a[i+1][i] = 0.0;
            for (j = i+1; j <= m; j++) {
                tmp       = cos_theta * a[i][j] - sin_theta * a[i+1][j];
                a[i+1][j] = sin_theta * a[i][j] + cos_theta * a[i+1][j];
                a[i][j]   = tmp;
            }

            // Q = P_i * Q^t----------------------------------------------
            for (j = 0; j <= m; j++) {
                tmp        = cos_theta * Qt[i][j] - sin_theta * Qt[i+1][j];
                Qt[i+1][j] = sin_theta * Qt[i][j] + cos_theta * Qt[i+1][j];
                Qt[i][j]   = tmp;
            }
        }

        // RQ���v�Z-------------------------------------------------------------
        for (i = 0; i <= m; i++) {
            for (j = 0; j <= m; j++)
                work[j] = a[i][j];
            for (j = 0; j <= m; j++) {
                tmp = 0.0;
                for (k = 0; k <= m; k++) {
                    tmp += work[k] * Qt[j][k];
                }
                a[i][j] = tmp;
            }
        }

        iter++;
    }
}

// ���C���֐�-----------------------------------------------------------------------
int main(void) {

    double a[N][N] = {
        { 16.0,  -1.0,   1.0,   2.0 },
        {  2.0,  12.0,   1.0,  -1.0 },
        {  0.0,   3.0, -24.0,   2.0 },
        {  0.0,   0.0,   1.0,  20.0 }
    };
    int i;

    qr_method(a);

    // ���ʏo��-----------------------------------------------------------------
    printf("�ŗL�l�͈ȉ��̒ʂ�:\n");
    for (i = 0; i < N; i++) {
        printf("%11.6f\n", a[i][i]);
    }

    return 0;
}
