#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#define N  3            // ������
#define EPS 1e-8        // ��������p
#define MAX_ITER 100    // �ő唽����

// ���������N���X�i���ϐ��Łj------------------------------------------------------
class BU_N {
public:

	// �f�[�^�����o---------------------------------------------------
	double v;
	double d[N] = { 0.0 };

	// �R���X�g���N�^-------------------------------------------------
	BU_N(double x = 0.0) : v(x) { }

	// �o�͊֐�-------------------------------------------------------
	void print() {
		printf("%lf < ", v);
		for (int i = 0; i < N - 1; i++) {
			printf("%lf, ", d[i]);
		}
		printf("%lf >\n", d[N - 1]);
	}
};

// ���������N���X�ɑ΂��鉉�Z�q���d��`--------------------------------------------
/* �P���}�C�i�X���Z�q */
BU_N operator-(const BU_N& a) {
	BU_N r(-a.v);
	for (int i = 0; i < N; i++) { r.d[i] = -a.d[i]; }
	return r;
}
/* ���Z */
BU_N operator+(const BU_N& a, const BU_N& b) {
	BU_N r(a.v + b.v);
	for (int i = 0; i < N; i++) { r.d[i] = a.d[i] + b.d[i]; }
	return r;
}
/* ���Z */
BU_N operator-(BU_N a, BU_N b) {
	BU_N r(a.v - b.v);
	for (int i = 0; i < N; i++) { r.d[i] = a.d[i] - b.d[i]; }
	return r;
}
/* ��Z */
BU_N operator*(const BU_N& a, const BU_N& b) {
	BU_N r(a.v * b.v);
	for (int i = 0; i < N; i++) { r.d[i] = a.d[i] * b.v + a.v * b.d[i]; }
	return r;
}
/* ���Z */
BU_N operator/(const BU_N& a, const BU_N& b) {
	BU_N r(a.v / b.v);
	for (int i = 0; i < N; i++) {
		r.d[i] = (a.d[i] * b.v - a.v * b.d[i]) / (b.v * b.v);
	}
	return r;
}

// �����֐�--------------------------------------------------------------------
BU_N sin(const BU_N& x) {
	BU_N r(std::sin(x.v));
	for (int i = 0; i < N; i++) { r.d[i] = std::cos(x.v) * x.d[i]; }
	return r;
}

// �����s�{�b�g�I��t���K�E�X�̏����@----------------------------------------------
void piboted_gauss(double a[][N], double b[], double x[])
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
		} else if (max_i != i) {

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

/* �ő�l�m���� */
double vector_norm_max(double x[]) {
	double max;
	int i;

	max = 0.0;
	for (i = 0; i < N; i++) {
		if (fabs(x[i]) > max) {
			max = fabs(x[i]);
		}
	}

	return max;
}

// BU_N�^�̔z���������------------------------------------------------------------
void init(BU_N BU_x[], const double x[]) {

	for (int i = 0; i < N; i++) {

		/* �x�N�g��x�̒l�����������N���X�̒l�ɑ�� */
		BU_x[i].v = x[i];

		/* ���R�r�s��̒l��P�ʍs��ŏ����� */
		for (int j = 0; j < N; j++) {
			BU_x[i].d[j] = (i == j) ? 1.0 : 0.0;
		}
	}
}

// BU_N�^�̔z�񂩂�֐��l�ƃ��R�r�s��𒊏o----------------------------------------
void split(const BU_N BU_x[], double x[], double d[][N]) {

	for (int i = 0; i < N; i++) {

		x[i] = BU_x[i].v;
		for (int j = 0; j < N; j++) {
			d[i][j] = BU_x[i].d[j];
		}
	}
}

// �֐��̒�`----------------------------------------------------------------------
void func(BU_N f[], BU_N x[]) {
	f[0] = 5.0 * x[0] * x[0] * x[0] - 2.0 * x[0] * x[0] + x[1] * x[2] - x[1] - 3.0;
	f[1] = 2.0 * x[1] * x[1] * x[2] + 7.0 * x[0] - 2.0 * x[2] - 8.0;
	f[2] = x[0] * x[1] * x[2] + 5.0 * x[0] * x[1] - 4.0 * x[2] + 1.0;
}

// �j���[�g���@�i���ϐ��Łj--------------------------------------------------------
void newton_n(void (*func)(BU_N f[], BU_N x[]), double x[], int* r) {

	BU_N BU_x[N], BU_y[N];
	double f[N], J[N][N], delta_x[N];

	int iter = 0;
	do {
		// �x�N�g��x�̒l�����������^�ϐ��֑��
		init(BU_x, x);

		// ���������^�Ŋ֐��l�ƕΔ����l���v�Z
		func(BU_y, BU_x);

		// �֐��l�ƕΔ����l�𕪗����Ē��o
		split(BU_y, f, J);

		// �ߎ����̏C���ʁ�x���v�Z
		piboted_gauss(J, f, delta_x);

		// �ߎ���x�̍X�V
		for (int i = 0; i < N; i++) x[i] -= delta_x[i];

		iter++;
	} while (vector_norm_max(delta_x) > EPS && iter < MAX_ITER);

	// �߂�l���Z�b�g
	*r = iter;
}

// ���C���֐�----------------------------------------------------------------------
int main() {

	double x[N] = { 1.0,1.0,1.0 };
	int iter;

	newton_n(func, x, &iter);

	// ���ʏo��--------------------------------------------------------
	if (iter == MAX_ITER) {
		printf("����������܂���ł����D\n");
		return 1;
	} else {
		printf("����������܂����D\n");
		for (int i = 0; i < N; i++) {
			printf("x[%d] = %lf\n", i, x[i]);
		}
		printf("�����񐔂�%d��ł��D\n", iter);
	}

	return 0;
}