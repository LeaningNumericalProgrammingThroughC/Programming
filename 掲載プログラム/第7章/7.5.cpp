#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#define EPS 1e-8       // ��������p臒l
#define MAX_ITER 7   // �ő唽����

// ���f���N���X----------------------------------------------------------------
class Complex {

	// �f�[�^�����o------------------------------------------------------
	double Re, Im;

public:
	/* �R���X�g���N�^ */
	Complex(double x = 0, double y = 0) : Re(x), Im(y) { }

	/* �o�͊֐� */
	void print() {
		printf("%.6lf + %.6lfi ", Re, Im);
	}

	/* �f�[�^�ւ̃A�N�Z�X */
	double get_Re() const { return Re; }
	double get_Im() const { return Im; }

	/* ��Βl�̌v�Z */
	double abs() const {
		return sqrt(Re * Re + Im * Im);
	}

	// ���Z�q�I�[�o�[���[�h---------------------------------------------
	/* �}�C�i�X���Z�q */
	Complex operator-() const {
		return Complex(-Re, -Im);
	}
	/* ���Z */
	friend Complex operator+(const Complex& x, const Complex& y) {
		return Complex(x.Re + y.Re, x.Im + y.Im);
	}
	/* ���Z */
	friend Complex operator-(const Complex& x, const Complex& y) {
		return Complex(x.Re - y.Re, x.Im - y.Im);
	}
	/* ��Z */
	friend Complex operator*(const Complex& x, const Complex& y) {
		return Complex(x.Re * y.Re - x.Im * y.Im,
			x.Re * y.Im + x.Im * y.Re);
	}
	/* ���Z */
	friend Complex operator/(const Complex& x, const Complex& y) {
		double work = y.Re * y.Re + y.Im * y.Im;

		return Complex((x.Re * y.Re + x.Im * y.Im) / work,
			(-x.Re * y.Im + x.Im * y.Re) / work);
	}
	/* �w���֐� */
	friend Complex exp(const Complex& x) {
		double eR = std::exp(x.Re);
		return Complex(eR * cos(x.Im), eR * sin(x.Im));
	}
};

// ���������N���X�i���f���Łj---------------------------------------------
class BU_C {

	// �f�[�^�����o--------------------------------------------------
	Complex v;
	Complex d;

public:

	// �R���X�g���N�^------------------------------------------------
	BU_C(double v_real = 0.0, double v_image = 0.0,
		double d_real = 0.0, double d_image = 0.0) {
		v = Complex(v_real, v_image);
		d = Complex(d_real, d_image);
	}

	BU_C(Complex x, Complex y = 0.0) {
		v = x;
		d = y;
	}

	/* �֐��l�C�����l�̎擾 */
	Complex value() const { return v; }
	Complex derivative() const { return d; }

	/* �����l�̏����� */
	void init() {
		d = Complex(1.0, 0.0);
	}

	/* �o�͊֐� */
	void print() {
		v.print(); printf("< "); d.print(); printf(">");
	}

	// ���Z�q�I�[�o�[���[�h------------------------------------------
	/* �P���}�C�i�X���Z�q */
	BU_C operator-() const { return BU_C(-v, -d); }

	/* ���Z */
	friend BU_C operator+(const BU_C& x, const BU_C& y) {
		return BU_C(x.v + y.v, x.d + y.d);
	}
	/* ���Z */
	friend BU_C operator-(const BU_C& x, const BU_C& y) {
		return BU_C(x.v - y.v, x.d - y.d);
	}
	/* ��Z */
	friend BU_C operator*(const BU_C& x, const BU_C& y) {
		return BU_C(x.v * y.v, x.d * y.v + x.v * y.d);
	}
	/* ���Z */
	friend BU_C operator/(const BU_C& x, const BU_C& y) {
		return BU_C(x.v / y.v,
			(x.d * y.v - x.v * y.d) / (y.v * y.v));
	}

	// �����֐�------------------------------------------------------
	friend BU_C exp(const BU_C& x) {
		Complex e = exp(x.v);
		return BU_C(e, e * x.d);
	}
};

// �֐��̒�`--------------------------------------------------------------
BU_C func(BU_C x) {

	return 3.0*x*x*x + x*x + 2.0;
}

// �j���[�g���@-----------------------------------------------------------------
Complex newton_C(BU_C(*f)(BU_C), Complex x_0, int* result, int* iterations) {

	BU_C x(x_0), y;
	Complex delta_x;
	int iter = 0, converged = 0;

	while (iter < MAX_ITER) {
		x.init();

		y = f(x);
		delta_x = y.value() / y.derivative();
		x = x - delta_x;

		iter++;

		// ��������--------------------------------------------------
		if (delta_x.abs() <= EPS) {
			converged = 1;
			break;
		}
	}

	// �߂�l�̃Z�b�g---------------------------------------------------
	*iterations = iter;
	*result = converged;
	return x.value();
}

// ���C���֐�-------------------------------------------------------------------
int main() {

	Complex y;
	int iter, result;

	y = newton_C(func, Complex(1.0, 1.0), &result, &iter);

	// ���ʕ\��--------------------------------------------------------
	if (result == 1) {
		printf("�����񐔂�%d��ł��D\n", iter);
		printf("���Fx = "); y.print(); printf("\n");
	} else {
		printf("�������܂���ł����D\n");
	}

	return 0;
}