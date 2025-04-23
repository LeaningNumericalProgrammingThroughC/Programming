#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>

#define DIM_Y 2                       // y�̎�����
#define DIM_Z (DIM_Y + DIM_Y*DIM_Y)   // ���R�r�s����܂߂�������
#define N_DIV 100                     // �����Q�E�N�b�^�@�̕�����
#define EPS 1.0e-8                    // �j���[�g���@�̎�������p臒l
#define MAX_ITER  50                  // �j���[�g���@�̍ő唽����

#define N_DIM 6        // ������

// �ő�l�m����------------------------------------------------
double vector_norm_max(double x[]) {

	double max;
	int i;

	max = fabs(x[0]);
	for (i = 1; i < DIM_Y; i++) {
		if (fabs(x[i]) > max) {
			max = fabs(x[i]);
		}
	}

	return max;
}

// �����s�{�b�g�I��t���K�E�X�̏����@----------------------------------------------
void pivoted_gauss(double a[][DIM_Y], double b[], double x[])
{
	int i, j, k, max_i;
	double m, s, max, tmp;

	// �O�i����------------------------------------------------------------
	for (i = 0; i < DIM_Y - 1; i++)
	{
		// �����s�{�b�g�I��-----------------------------------------
		max = fabs(a[i][i]);
		max_i = i;
		for (j = i + 1; j < DIM_Y; j++) {

			if (fabs(a[j][i]) > max) {
				max = fabs(a[j][i]);
				max_i = j;
			}
		}

		if (max < EPS) {

			printf("�W���s�񂪐����ł͂���܂���");
			exit(1);
		} else if (max_i != i) {

			for (j = i; j < DIM_Y; j++) {
				tmp = a[i][j]; a[i][j] = a[max_i][j]; a[max_i][j] = tmp;
			}
			tmp = b[i]; b[i] = b[max_i]; b[max_i] = tmp;
		}

		// �O�i����-------------------------------------------------
		for (j = i + 1; j < DIM_Y; j++) {

			m = a[j][i] / a[i][i];

			for (k = i + 1; k < DIM_Y; k++) {
				a[j][k] = a[j][k] - m * a[i][k];
			}
			b[j] = b[j] - m * b[i];
		}
	}

	// ��ޑ��------------------------------------------------------------
	x[DIM_Y - 1] = b[DIM_Y - 1] / a[DIM_Y - 1][DIM_Y - 1];
	for (i = DIM_Y - 2; i >= 0; i--) {

		s = 0;
		for (j = i + 1; j < DIM_Y; j++) {
			s += a[i][j] * x[j];
		}
		x[i] = (b[i] - s) / a[i][i];
	}
}

// ���ϐ������Q�E�N�b�^�@------------------------------------------------
void runge_kutta_system(double a, double b, int n, double y[],
						       void (*f)(double, double[], double[])) {
	double x = a, h = (b - a) / n;
	double k1[N_DIM], k2[N_DIM], k3[N_DIM], k4[N_DIM], tmp[N_DIM];
	int i, j;

	// ���o����1�s�ڂ��o��-----------------------------------------
	/*
	printf("    x         ");
	for (j = 0; j < N_DIM; j++) printf("y[%d]        ", j);
	printf("\n");

	printf("%7.3lf   ", x);
	for (j = 0; j < N_DIM; j++) printf("%10.5lf  ", y[j]);
	printf("\n");
	*/

	for (i = 0; i < n; i++) {

		// k1 ���v�Z -----------------------------------------
		f(x, y, k1);

		// k2 ���v�Z -----------------------------------------
		for (j = 0; j < N_DIM; j++) tmp[j] = y[j] + h * k1[j] / 2.0;
		f(x + h / 2.0, tmp, k2);

		// k3 ���v�Z -----------------------------------------
		for (j = 0; j < N_DIM; j++) tmp[j] = y[j] + h * k2[j] / 2.0;
		f(x + h / 2.0, tmp, k3);

		// k4 ���v�Z -----------------------------------------
		for (j = 0; j < N_DIM; j++) tmp[j] = y[j] + h * k3[j];
		f(x + h, tmp, k4);

		// y ���v�Z	------------------------------------------
		for (j = 0; j < N_DIM; j++)
			y[j] = y[j] + h / 6.0
			* (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);

		x += h;

		// ���ʏo��	------------------------------------------
		/*
		printf("%7.3lf   ", x);
		for (j = 0; j < N_DIM; j++) printf("%10.5lf  ", y[j]);
		printf("\n");
		*/
	}
}

// �����������̉E�ӂ��`---------------------------------------------------------
void f(double x, double y[], double r[]) {

	r[0] = y[1];
	r[1] = -0.1 * y[1] - y[0] * y[0] * y[0] + 0.3 * cos(x);
	r[2] = y[3];
	r[3] = -3.0 * y[0] * y[0] * y[2] - 0.1 * y[3];
	r[4] = y[5];
	r[5] = -3.0 * y[0] * y[0] * y[4] - 0.1 * y[5];
}

// �V���[�e�B���O�@---------------------------------------------------------------
void shooting(double a, double b, double y_a[], 
	           void (*func)(double, double[], double[]), int* r) {

	double F[DIM_Y], dF[DIM_Y][DIM_Y], dy[DIM_Y], z[DIM_Z];
	int iter = 0;

	// �j���[�g���@�����s--------------------------------------------------
	while (iter++ < MAX_ITER ) {

		// �e�֐��̏����l���Z�b�g
		z[0] = y_a[0];
		z[1] = y_a[1];
		z[2] = 1; z[3] = 0;
		z[4] = 0; z[5] = 1;

		// �����Q�E�N�b�^�@��y(b)�ƃ��R�r�s��d��(b,��)/d�����v�Z
		runge_kutta_system(a, b, N_DIV, z, func);

		// F(��)���v�Z
		F[0] = z[0] - y_a[0];
		F[1] = z[1] - y_a[1];

		// ���R�r�s��dF/d�����v�Z
		dF[0][0] = z[2] - 1.0;
		dF[0][1] = z[4];
		dF[1][0] = z[3];
		dF[1][1] = z[5] - 1.0;

		// �j���[�g���@�̌����ɏ]��F(��)�̕������t�ɂ���
		for (int j = 0; j < DIM_Y; j++) F[j] = -F[j];

		// �K�E�X�̏����@�Ńj���[�g���@�̏C���ʂ����߂�
		pivoted_gauss(dF, F, dy);

		// �ߎ����̍X�V
		for (int j = 0; j < DIM_Y; j++) y_a[j] += dy[j];

		// ��������
		if (vector_norm_max(dy) <= EPS) break;
	}

	// �߂�l���Z�b�g
	*r = iter;
}

// ���C���֐�---------------------------------------------------------------------
int main() {

	double a = 0.0, b = 2.0 * M_PI, y_a[DIM_Y];
	int iter = 0;

	// �����l�̐ݒ�
	y_a[0] = 1.0;
	y_a[1] = 1.0;

	shooting(a, b, y_a, f, &iter);

	// ���ʏo��------------------------------------------------------------
	if (iter > MAX_ITER) {
		printf("�����񐔂̏���ɒB���܂����D\n");
	} else {
		printf("���̏����l�͈ȉ��̂Ƃ���D�i�����񐔁F%d��j\n", iter);
		for (int j = 0; j < DIM_Y; j++) {
			printf("y_a[%d] = %8.6lf\n", j, y_a[j]);
		}
	}

	return 0;
}