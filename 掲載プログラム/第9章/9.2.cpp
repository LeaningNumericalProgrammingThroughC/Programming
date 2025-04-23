#include <stdio.h>
#include <math.h>

#define N    4            // ������
#define EPS  1e-8         // ��������p臒l
#define MAX_ITER  100     // �ő唽����

// ���R�r�@-----------------------------------------------------------------------------
void jacobi_method(double A[][N], double U[][N], int* r){

    double max_val, A_p_p, A_p_q, A_q_q;
    double sin_2_theta, cos_2_theta, tan_2_theta, sin_theta, cos_theta, tmp;
    int i, j, p, q, iter;
    
    // U ��P�ʍs��ɏ�����---------------------------------------------------------
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i == j) {
                U[i][j] = 1.0;
            } else {
                U[i][j] = 0.0;
            }
        }
    }

    // ��������---------------------------------------------------------------------
    for (iter = 0; iter < MAX_ITER; iter++) {

        // ��Βl�ő�̔�Ίp�v�f��T��
        max_val = 0.0;
        p = 0;
        q = 1;
        for (i = 0; i < N; i++) {
            for (j = i + 1; j < N; j++) {
                if (fabs(A[i][j]) > max_val) {
                    max_val = fabs(A[i][j]);
                    p = i;
                    q = j;
                }
            }
        }

        // ��������F��Ίp�������\���������Ȃ�I��-------------------------------
        if (max_val < EPS) {
            break;
        }

        // �s��̐������㏑�������O�ɁC�v�Z�̌��ƂȂ鐬����ޔ�-----------------
        A_p_p = A[p][p];
        A_p_q = A[p][q];
        A_q_q = A[q][q];

        // sin�Ƃ�cos�Ƃ̌v�Z-----------------------------------------------------
        tan_2_theta = -2 * A_p_q / (A_p_p - A_q_q);
        tmp = sqrt(1 + tan_2_theta * tan_2_theta);
        cos_2_theta = 1.0 / tmp;
        sin_2_theta = tan_2_theta / tmp;
        cos_theta = sqrt(0.5 * (1.0 + cos_2_theta));
        sin_theta = sin_2_theta / (2.0 * cos_theta);

        // �s��A��p�s��q�s�̒l���X�V----------------------------------------------
        for (i = 0; i < N; i++) {
            tmp = A[q][i];
            A[q][i] = A[p][i] * sin_theta + tmp * cos_theta;
            A[p][i] = A[p][i] * cos_theta - tmp * sin_theta;
        }

        // �s��A��p�s��q�s�̒l��p���q��ɃR�s�[----------------------------------
        for (i = 0; i < N; i++) {
            A[i][p] = A[p][i];
            A[i][q] = A[q][i];
        }

        // A_pp,A_pq,A_qp,A_qq�����̍X�V------------------------------------------
        tmp = 2.0 * A_p_q * sin_theta * cos_theta;
        A[p][p] = A_p_p *cos_theta *cos_theta + A_q_q *sin_theta *sin_theta - tmp;
        A[q][q] = A_p_p *sin_theta *sin_theta + A_q_q *cos_theta *cos_theta + tmp;
        A[p][q] = A[q][p] = 0;

        // �ŗL�x�N�g���s�� U �̍X�V----------------------------------------------
        for (i = 0; i < N; i++) {
            tmp = U[i][p];
            U[i][p] = cos_theta * tmp - sin_theta * U[i][q];
            U[i][q] = sin_theta * tmp + cos_theta * U[i][q];
        }
    }

    *r = iter;
}

// ���C���֐�-----------------------------------------------------------------------------
int main(void){

    double A[N][N] = {
        {  5.0,  1.0, -3.0, -4.0 },
        {  1.0,  2.0, -5.0,  2.0 },
        { -3.0, -5.0,  7.0,  3.0 },
        { -4.0,  2.0,  3.0, -6.0 }
    };
    double U[N][N]; 
    int i, j, iter;

    jacobi_method(A, U, &iter);

    // ���ʂ̕\��-----------------------------------------------------------------
    if (iter == MAX_ITER) {
        printf("�����񐔂�����ɒB���܂����D\n");
    } else {
        printf("�ŗL�l:\n");
        for (i = 0; i < N; i++) {
            printf("%10.6f ", A[i][i]);
        }
        printf("\n\n�ŗL�x�N�g�� (U �̊e��):\n");
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                printf("%10.6f ", U[i][j]);
            }
            printf("\n");
        }
    }

    return 0;
}
