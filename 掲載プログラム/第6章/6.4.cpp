#include <stdio.h>
#include <cmath>

// ���������N���X�i�{�g���A�b�v�^�j--------------------------------------------
class BU {

	// �f�[�^�����o---------------------------------------------------
    double v;
    double d;

public:
	/* �R���X�g���N�^ */
	BU(double a = 0.0, double b = 0.0) : v(a), d(b) { }

    /* �֐��l�C�����l�̎擾 */
	double value() const { return v; }
	double derivative() const { return d; }

	/* �����l�̏����� */
	void init() {
		d = 1.0;
	}

	/* �o�͊֐��C �l < �����l > �̌`���ŕ\�� */
	void print() {
		printf("%lf < %lf >", v, d);
	}

	// ���Z�q�I�[�o�[���[�h-------------------------------------------
	/* �P���}�C�i�X���Z�q */
	BU operator-() const { return BU(-v, -d); }

	/* ���Z */
	friend BU operator+(const BU& x, const BU& y) {
		return BU(x.v + y.v, x.d + y.d);
	}
	/* ���Z */
	friend BU operator-(const BU& x, const BU& y) {
		return BU(x.v - y.v, x.d - y.d);
	}
	/* ��Z */
	friend BU operator*(const BU& x, const BU& y) {
		return BU(x.v * y.v, x.d * y.v + x.v * y.d);
	}
	/* ���Z */
	friend BU operator/(const BU& x, const BU& y) {
		return BU( x.v / y.v, (x.d*y.v  - x.v*y.d) / (y.v*y.v) );
	}

	// �����֐�-------------------------------------------------------
	friend BU sqrt(const BU& x){
		double s = std::sqrt(x.v);
		return BU(s, 0.5 * x.d / s);
	}
	friend BU sin(const BU& x) {
		return BU(std::sin(x.v), std::cos(x.v) * x.d);
	}
	friend BU cos(const BU& x) {
		return BU(std::cos(x.v), -std::sin(x.v) * x.d);
	}
	friend BU exp(const BU& x){
		double e = std::exp(x.v);
		return BU(e, e * x.d);
	}
	friend BU log(const BU& x){
		return BU(std::log(x.v), x.d / x.v);
	}
};

// ���C���֐�------------------------------------------------------------------
int main() {

	BU x(2.0,1.0), y;

	y =  sin(x * x * x - exp(x)) / (x * x + 2.0);

	printf("y = ");	y.print(); printf("\n");

	return 0;
}