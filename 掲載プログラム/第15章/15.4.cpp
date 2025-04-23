#include <stdio.h>
#define _USE_MATH_DEFINES 
#include <cmath>

#define N 2            // ������
#define DIV 100        // ������
#define EPS 1.0e-8     // ��������p臒l
#define MAX_ITER  50   // �ő唽����

// �ő�l�m����------------------------------------------------
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

// �����������̉E��-------------------------------------------------------------
void func(BU_N f[], double x, BU_N y[]) {

	f[0] = y[1];
	f[1] = -0.1 * y[1] - y[0] * y[0] * y[0] + 0.3 * cos(x);
}

// ���������^�p�̃����Q�E�N�b�^�@-----------------------------------------------
void runge_kutta_system_BU(BU_N f_end[N], BU_N y_0[], double start, double end, 
	void (*func)(BU_N[], double, BU_N[])) {

	BU_N y[N], y_old[N], f[N], y_tmp[N];
	BU_N k1[N], k2[N], k3[N], k4[N];
	double h = (end - start) / DIV, x = start;

	int k = 0;

	for (int i = 0; i < N; i++) {
		y[i] = y_0[i];
	}

	for (k = 0; k < DIV; k++) {

		for (int j = 0; j < N; j++)	y_old[j] = y[j];

		func(f, x, y);
		for (int j = 0; j < N; j++)  k1[j] = f[j];

		for (int j = 0; j < N; j++) y_tmp[j] = y[j] + h / 2. * k1[j];
		func(f, x + h / 2, y_tmp);
		for (int j = 0; j < N; j++) k2[j] = f[j];

		for (int j = 0; j < N; j++) y_tmp[j] = y[j] + h / 2. * k2[j];
		func(f, x + h / 2, y_tmp);
		for (int j = 0; j < N; j++) k3[j] = f[j];

		for (int j = 0; j < N; j++) y_tmp[j] = y[j] + h * k3[j];
		func(f, x + h, y_tmp);
		for (int j = 0; j < N; j++) k4[j] = f[j];

		for (int j = 0; j < N; j++) {
			y[j] = y[j] + (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
		}

		x += h;
	}

	for (int j = 0; j < N; j++) f_end[j] = y[j];
}

// �V���[�e�B���O�@---------------------------------------------------------------
int shooting_BU(double alpha[], double start, double end,
	void (*func)(BU_N[], double, BU_N[]), int* r) {

	BU_N ad_g[N], ad_alpha[N];
	double g[N], dg[N][N], d_alpha[N];
	int converged = 0, n_iter = 0, i, j;


	while (n_iter < MAX_ITER) {

		// alpha �����������^�ɕϊ�----------------------------------------------
		init(ad_alpha, alpha);

		// �����Q�E�N�b�^�@��G(��)���v�Z---------------------------------------------
		runge_kutta_system_BU(ad_g, ad_alpha, start, end, func);
		for (i = 0; i < N; i++) ad_g[i] = ad_g[i] - ad_alpha[i];

		// ���������^ad_g���֐��lg(��)�ƃ��R�r�s��dg�ɕ���-------------------------
		split(ad_g, g, dg);

		// �j���[�g���@�̌����ɏ]��g(��)�̕������t�ɂ���------------------------
		for (j = 0; j < N; j++) g[j] = -g[j];

		// �K�E�X�̏����@�ɂ��C�j���[�g���@�̏C���ʂ����߂�--------------------
		pivoted_gauss(dg, g, d_alpha);

		// �ߎ����̍X�V---------------------------------------------------------
		for (int j = 0; j < N; j++) alpha[j] += d_alpha[j];

		n_iter++;

		if (vector_norm_max(d_alpha) <= EPS) {
			converged = 1;
			break;
		}
	}

	// �߂�l���Z�b�g
	*r = n_iter;
	return converged;
}

// ���C���֐�---------------------------------------------------------------------
int main() {

	double a = 0.0, b = 2.0 * M_PI, y_a[N];
	int converged, iter = 0;

	// �����l�̐ݒ�
	y_a[0] = 1.0;
	y_a[1] = 1.0;

	converged = shooting_BU(y_a, a, b, func, &iter);
	
	// ���ʏo��------------------------------------------------------------
	if (converged) {
		printf("���̏����l�͈ȉ��̂Ƃ���D�i�����񐔁F%d��j\n", iter);
		for (int j = 0; j < N; j++) {
			printf("y_a[%d] = %8.6lf\n", j, y_a[j]);
		}
	} else {
		printf("�����񐔂̏���ɒB���܂����D\n");
	}

	return 0;
}
