#include <stdio.h>
#include <cmath>

#define N 3    // ������

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
	for (int i = 0; i < N; i++) { r.d[i] = a.d[i]*b.v + a.v*b.d[i]; }
	return r;
}
/* ���Z */
BU_N operator/(const BU_N& a, const BU_N& b) {
	BU_N r(a.v / b.v);
	for (int i = 0; i < N; i++) { 
		r.d[i] = ( a.d[i]*b.v - a.v*b.d[i] ) / (b.v * b.v); 
	}
	return r;
}

// �����֐�--------------------------------------------------------------------
BU_N sin(const BU_N& x) {
	BU_N r(std::sin(x.v));
	for (int i = 0; i < N; i++) { r.d[i] = std::cos(x.v) * x.d[i]; }
	return r;
}

   // �E�E�E�@�i�ȉ��Ccos(x),exp(x)�������l�j


// ���C���֐�------------------------------------------------------------------
int main() {

	double w[N] = { 1,0,0 };

	BU_N x(2.0),y(3.0),z(4.0);

	/* �e�ϐ��̕Δ����l�̏����� */
	x.d[0] = y.d[1] = z.d[2] = 1.0;

	BU_N f =  5.0 * x * x - 2.0 * x * y * y + 3.0 * y * z;
	
	printf("f <df/dx, df/dy, df/dz> = ");	f.print(); printf("\n");

	return 0;
}
