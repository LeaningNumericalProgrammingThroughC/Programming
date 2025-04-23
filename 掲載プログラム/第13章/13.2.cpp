#include <stdio.h>
#define _USE_MATH_DEFINES  // Visual Studio����M_PI���g���ꍇ�K�v
#include <math.h>

#define N 100    // �W�{�_��

// �֐��̒�`----------------------------------------------------------------
double func(double x) {
    return x + cos(x);
}

// �K�E�X�E�`�F�r�V�F�t�ϕ�����----------------------------------------------
double gauss_chebyshev(double (*f)(double), int n) {

    double s = 0.0, theta;
    int i;

    for (i = 1; i <= n; i++) {
        theta = (2.0 * i - 1.0) * M_PI / (2.0 * n);
        s += f(cos(theta)) * sin(theta);
    }

    return M_PI / n * s;
}

// ���C���֐�----------------------------------------------------------------
int main() {

    printf("�K�E�X�E�`�F�r�V�F�t�ϕ�����: %12.10lf\n", gauss_chebyshev(func, N));

    return 0;
}