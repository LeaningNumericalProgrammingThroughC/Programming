#include <stdio.h>
#include <math.h>

#define EPS 1e-8         // ��������p臒l
#define MAX_ITER 50      // �ő唽����

// ���������N���X�i�{�g���A�b�v�^�j--------------------------------------------
class BU {

	double v;
	double d;

public:

	/* �R���X�g���N�^ */
	BU(double a = 0.0, double b = 0.0) {
		v = a; d = b;
	}
	/* �֐��l�C�����l�̎擾 */
	double value() const { return v; }
	double derivative() const { return d; }

	/* �����l�̏����� */
	void init() {
		d = 1.0;
	}
	/* �o�͊֐� */
	void print() {
		printf("%lf < %lf >", v, d);
	}

	/* �}�C�i�X���Z�q */
	BU operator-() const { return BU(-v, -d); }

	/* ���Z (BU �{ BU) */
	friend BU operator+(const BU& x, const BU& y) {
		return BU(x.v + y.v, x.d + y.d);
	}
	/* ���Z (BU - BU) */
	friend BU operator-(const BU& x, const BU& y) {
		return BU(x.v - y.v, x.d - y.d);
	}
	/* ��Z (BU * BU) */
	friend BU operator*(const BU& x, const BU& y) {
		return BU(x.v * y.v, x.d * y.v + x.v * y.d);
	}
	/* ���Z (BU / BU) */
	friend BU operator/(const BU& x, const BU& y) {
		return BU(x.v / y.v, (x.d * y.v - x.v * y.d) / (y.v * y.v));
	}

	/* �����֐� */
	friend BU sqrt(const BU& x) {
		return BU(sqrt(x.v), 0.5 * x.d / (sqrt(x.v)));
	}
	friend BU sin(const BU& x) {
		return BU(sin(x.v), cos(x.v) * x.d);
	}
	friend BU cos(const BU& x) {
		return BU(cos(x.v), -sin(x.v) * x.d);
	}
	friend BU exp(const BU& x) {
		return BU(exp(x.v), exp(x.v) * x.d);
	}
	friend BU log(const BU& x) {
		return BU(log(x.v), x.d / x.v);
	}

};

// �֐��̒�`-------------------------------------------------------------------
BU func(BU x) {

	return x * x * x - 3 * x * x - exp(2 * x - 5) + 3;
}

// �j���[�g���@-----------------------------------------------------------------
double newton( BU (*f)(BU), double x_0, int *result, int* iterations) {

	BU x(x_0), y;
	double delta_x;
	int iter = 0, converged = 0;

	while (iter < MAX_ITER) {

		// �����l���������i0���Z�b�g�j
		x.init();  

		y = f(x);
		delta_x = y.value() / y.derivative();
		x = x - delta_x;
		
		iter++;

		// ��������
		if (fabs(delta_x) <= EPS) {
			converged = 1; 
			break;
		}
	}

	// �߂�l�̃Z�b�g
	*iterations = iter;
	*result = converged;
	return x.value();
}

// ���C���֐�-------------------------------------------------------------------
int main() {

	double y;
	int iter, result;
	
	y = newton(func, 1.0, &result, &iter);
	
	// ���ʕ\��--------------------------------------------------------
	if (result == 1) {
		printf("�����񐔂�%d��C����%lf�ł��D\n", iter, y);
	} else {
		printf("�������܂���ł����D\n");
	}

	return 0;
}