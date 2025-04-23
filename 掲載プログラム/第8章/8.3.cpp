#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#define N    3            // �������ix:2���� + t:1�����j
#define DELTA_S    0.1    // ���Ȑ��ǐՂ̍��ݕ�
#define EPS   1e-10       // �j���[�g���@�̎�������p
#define MAX_ITER 100      // �j���[�g���@�̍ő唽����
#define MAX_STEP 1000     // �z���g�s�[�@�̍ő�X�e�b�v��

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

// �O���[�o���ϐ�--------------------------------------------------------------------

// g_f_gamma�Ff(��)
double g_f_gamma_1, g_f_gamma_2;
// g_n�F���Ȑ��̐ڃx�N�g���C������q(x,t)=0�̖@���x�N�g��
double g_n[N];
// g_p_c�F�}8.4�̓_p_c�C�}8.4�̓_p_c�C�j���[�g���@�̏����_
double g_p_c[N];

// �֐�f(x)�̑�1�����idouble�p�j-----------------------------------------------------
double f_1_double(double x_1, double x_2) {
	return  2.0 * x_1 * x_1 * x_1 + 2. * x_1 * x_2 - 22.0 * x_1 + x_2 * x_2 + 13.0;
}
// �֐�f(x)�̑�2�����idouble�p�j-----------------------------------------------------
double f_2_double(double x_1, double x_2) {
	return  x_1 * x_1 + 2. * x_1 * x_2 + 2. * x_2 * x_2 * x_2 - 14. * x_2 + 9.;
}

// �֐�f(x)�̑�1�����i���������^�p�j-------------------------------------------------
BU_N f_1(BU_N x_1, BU_N x_2) {
	return  2.0 * x_1 * x_1 * x_1 + 2. * x_1 * x_2 - 22.0 * x_1 + x_2 * x_2 + 13.0;
}
// �֐�f(x)�̑�2�����i���������^�p�j-------------------------------------------------
BU_N f_2(BU_N x_1, BU_N x_2) {
	return  x_1 * x_1 + 2. * x_1 * x_2 + 2. * x_2 * x_2 * x_2 - 14. * x_2 + 9.;
}

// �z���g�s�[�֐�h(x)----------------------------------------------------------------
void func(BU_N h[], BU_N x[]) {

	BU_N BU_f[N - 1];

	// x[0]�͎�����x_1,  x[1]�Fx_2 , x[2]�F t

	// �z���g�s�[������
	h[0] = f_1(x[0], x[1]) + (x[2] - 1.0) * g_f_gamma_1;
	h[1] = f_2(x[0], x[1]) + (x[2] - 1.0) * g_f_gamma_2;

	// ����q(x)=0�̕�����
	h[2] = g_n[0] * (x[0] - g_p_c[0])
		+ g_n[1] * (x[1] - g_p_c[1])
		+ g_n[2] * (x[2] - g_p_c[2]);
}

// �j���[�g���@�i���ϐ��Łj----------------------------------------------------------
void newton_n(void (*func)(BU_N [], BU_N []), double x[], int* r) {

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

// �z���g�s�[�@----------------------------------------------------------------------
int homotopy(void (*func)(BU_N [], BU_N []), double x_s[]) {

	// f(��)�̒l���z���g�s�[�֐�h(x)�ɃZ�b�g����
	g_f_gamma_1 = f_1_double(x_s[0], x_s[1]);
	g_f_gamma_2 = f_2_double(x_s[0], x_s[1]);

	// ���Ȑ���̒ǐՓ_�i�X�V��j
	double p_new[N] = { x_s[0], x_s[1] , 0 };
	// ���Ȑ���̒ǐՓ_�i�X�V�O�j
	double p_old[N] = { x_s[0], x_s[1], -DELTA_S };

	int i, iter_n, iter_h = 0, num = 0;

	while (fabs( p_new[2] ) < 5.0 && iter_h < MAX_STEP ) {

		// ���Ȑ��̐ڐ������x�N�g�����v�Z
		for (i = 0; i < N; i++) {
			g_n[i] = p_new[i] - p_old[i]; 
		}

		// Newton�@�̏����_p_c���v�Z
		for (i = 0; i < N; i++) {
			g_p_c[i] = p_new[i] + DELTA_S * g_n[i] / vector_norm_max(g_n);
		}

		// ���̔����̏����Ƃ��čX�V��̒ǐՓ_���X�V�O�̒ǐՓ_�ɃZ�b�g
		for (i = 0; i < N; i++) p_old[i] = p_new[i];

		// p_new��Newton�@�̏����l���Z�b�g
		for (i = 0; i < N; i++) p_new[i] = g_p_c[i];

		// Newton�@�����s */
		newton_n(func, p_new, &iter_n);

		// �j���[�g���@�̔��U���`�F�b�N
		if (iter_n == MAX_ITER) {
			printf("���Ȑ��̒ǐՂɎ��s���܂����D\n");
			return 1;
		}

		// ���Ȑ�������t=0�ƌ������Ă���΋ߎ������o��
		if ((p_new[2] - 1.0) * (p_old[2] - 1.0) < 0 || p_new[2] == 1.0) {
			printf("�� ��_%d�F\n", ++num );
			for (i = 0; i < N -1; i++) {
				printf("x[%d]=%lf\n", i, ( fabs(p_new[2] -1.0) *p_old[i] 
					+ fabs(1.0 -p_old[2]) *p_new[i] ) / fabs(p_new[2] -p_old[2]) );
			}
			printf("\n");
		}

		iter_h++;
	}

	if (iter_h == MAX_STEP) {
		printf("�����񐔂�����ɒB���܂����D\n");
	} else {
		printf("t�̒l���ݒ�͈͂��z���܂����D\n");
	}

	return 0;
}

// ���C���֐�------------------------------------------------------------------------
int main() {

	double x_s[N - 1] = { -6.0, -6.0 };

	homotopy(func, x_s);

	return 0;
}