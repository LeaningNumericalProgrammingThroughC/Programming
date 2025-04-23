#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES  // Visual Studio���ŃR���p�C������ꍇ�K�v
#include <cmath>

#define N 5                    // ����
#define EPS_BISECTION 1e-1     // 2���@�p�̎�������臒l
#define EPS_DKA 1e-8           // DKA�@�p�̎�������臒l
#define MAX_ITER 100           // �ő唽����


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

// ������P(z)�̌v�Z----------------------------------------------------------
Complex P(double a[], Complex x) {
    int i;
    Complex r;

    r = a[0];
    for (i = 1; i <= N; i++) {
        r = r * x + a[i];
    }
    return r;
}

// DK�@�̏C����delta_x�̕���̌v�Z-------------------------------------------
Complex PI(Complex x[], int i) {
    int j;
    Complex r(1, 0);

    for (j = 0; j < N; j++) {
        if (i != j) {
            r = r * (x[i] - x[j]);
        }
    }
    return r;
}

// ������S�̌v�Z-------------------------------------------------------------
double S(double b[], double x) {
    int i;
    double r;

    r = b[0];
    for (i = 1; i <= N; i++) {
        r = r * x + b[i];
    }
    return r;
}

// 2���@---------------------------------------------------------------------
double bisection_DKA(double b[], double left, double right) {
    
    double center, S_center;
    double S_left = S(b, left);

    while ( (right - left) > EPS_BISECTION * right) {
        
        center = (right + left) / 2.0;
        S_center = S(b, center);

        if (S_center == 0.0) break;

        if (S_left * S_center < 0.0) {
            right = center;
        }
        else {
            left = center;
            S_left = S_center;
        }
    }

    return (left + right) / 2.0;
}

// DKA�@---------------------------------------------------------------------
int dka(double a[], Complex* z, int* r) {

    Complex delta_z;
    double R, b[N + 1], max_delta_z, zeta, eta = 0, tmp;
    int i, j, iter= 0, converged = 0;

    /* ���̏d�S���v�Z */
    zeta = -a[1] / (a[0] * N);

    for (i = 0; i < N + 1; i++)  b[i] = a[i];

     /* ������P^*�̌W�����v�Z */
    for (i = 0; i <= N-2; i++) {
        for (j = 1; j <= N - i; j++) {
            b[j] = b[j] + b[j - 1] * zeta;
        }
    }

    /* ������S�̌W�����v�Z */
    int m = 0;
    for (i = 2; i <= N; i++) {
        b[i] = -fabs(b[i]);
        if (fabs(b[i]) > 1.0e-8) m++;
    }

    /* �ł̌v�Z */
    for (i = 2; i <= N; i++) {
        tmp = pow(m * fabs(b[i]) / fabs(b[0]), 1.0 / (double)i);
        if (tmp > eta) eta = tmp;
    }

    /* 2���@�̎��s */
    R = bisection_DKA(b, 0, eta);
 
    /* �����l�̐ݒ� */
    for (j = 0; j < N; j++) {
        z[j] = zeta + R * exp(Complex(0, 1) * (2 * M_PI * ((double) j) / N)
            + M_PI / (2 * N));
    }

    /* DK�@�̔����X�V */
    while (iter < MAX_ITER) {
        max_delta_z = 0;
        for (i = 0; i < N; i++) {

            delta_z = P(a, z[i]) / (a[0] * PI(z, i));
            z[i] = z[i] - delta_z;

            if (delta_z.abs() > max_delta_z) {
                max_delta_z = delta_z.abs();
            }
        }
        iter++;

        // ��������---------------------------
        if (max_delta_z <= EPS_DKA) {
            converged = 1;
            break;
        }
    }

    /* �߂�l���Z�b�g */
    *r = iter;
    return converged;
}

// ���C���֐�---------------------------------------------------------------
int main(void) {

    Complex z[N];
    // ������P(z)�̌W���ia[i] �F x^(N-i)�̌W���j
    double a[N + 1] = { 1,-2,4,-5,3,1 };
    int i,iter, result;

    result = dka(a, z, &iter);

    // ���ʏo��-----------------------------------------------------
    if (result == 1) {
        printf("���͈ȉ��̒ʂ�ł��D�i�����񐔁F%d��j\n", iter);
        for (i = 0; i < N; i++) {
            printf("x[%d] = ", i); z[i].print(); printf("\n");
        }
    } else {
        printf("�����񐔂�����ɒB���܂����D\n");
    }

    return 0;
}


