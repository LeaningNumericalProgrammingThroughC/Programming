#include <stdio.h>
#include <math.h>

// �֐��̒�`-----------------------------------------------------------
double func(double x) {
    return x + cos(x);
}

// ��`����-------------------------------------------------------------
double trapezoidal(double (*f)(double), double a, double b, int n) {

    double s, h = (b - a) / n;
    int i;

    s = f(a) + f(b);
    for (i = 1; i < n; i++) {
        s += 2.0 * f(a + i * h);
    }

    return s * h / 2.0;
}

// �V���v�\������-------------------------------------------------------
double simpson(double (*f)(double), double a, double b, int n) {

    double s, h = (b - a) / (2.0 * n);
    int i;

    s = f(a) + 4.0*f(a + (2*n - 1)*h) + f(b);
    for (i = 1; i < n; i++) {
        s += 4.0 * f(a + (2*i-1)*h) + 2.0 * f(a + 2*i*h);
    }

    return s * h / 3.0;
}

// ���C���֐�-----------------------------------------------------------
int main() {

    double a = -1, b = 1;
    int n = 100;

    // ��`����
    printf("��`�����@�@�@ �F %12.10lf\n", trapezoidal(func, a, b, n));
    // �V���v�\������
    printf("�V���v�\������ �F %12.10lf\n", simpson(func, a, b, n));

    return 0;
}
