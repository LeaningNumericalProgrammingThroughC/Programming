#include <stdio.h>

#define N 10    // ������

// ���������� dy/dx = f(x, y) ���`--------------------------------
double f(double x, double y) {

    return 8 * x * y;
}

// �I�C���[�@-------------------------------------------------------
void euler(double a, double b, double y[], int n, 
                                  double (*f) (double, double)) {
    double x = a, h = (b - a) / n;

    for (int i = 0; i < n; i++) {
        y[i+1] = y[i] + h * f(x, y[i]);
        x = x + h;
    }
}

// �z�C���@---------------------------------------------------------
void heun(double a, double b, double y[], int n,
                                  double (*f) (double, double)) {
    double k1, k2, x = a, h = (b - a) / n;

    for (int i = 0; i < n; i++) {
        k1 = f(x, y[i]);
        k2 = f(x + h, y[i] + h * k1);

        y[i+1] = y[i] + h * (k1 + k2) / 2.0;
        x = x + h;
    }
}

// ���_�@-----------------------------------------------------------
void mid_point(double a, double b, double y[], int n,
                                  double (*f) (double, double)) {
    double k1, k2, x = a, h = (b - a) / n;

    for (int i = 0; i < n; i++) {
        k1 = f(x, y[i]);
        k2 = f(x + h / 2.0, y[i] + (h / 2.0) * k1);

        y[i+1] = y[i] + h * k2;
        x = x + h;
    }
}

// 3���̃����Q�E�N�b�^�@--------------------------------------------
void runge_kutta_3(double a, double b, double y[], int n,
                                  double (*f) (double, double)) {
    double k1, k2, k3, x = a, h = (b - a) / n;

    for (int i = 0; i < n; i++) {
        k1 = f(x, y[i]);
        k2 = f(x + h / 2.0, y[i] + h * k1 / 2.0);
        k3 = f(x + h, y[i] - h*k1 + 2.0*h*k2);

        y[i+1] = y[i] + h * (k1 + 4.0*k2 + k3) / 6.0;
        x = x + h;
    }
}

// 4���̃����Q�E�N�b�^�@--------------------------------------------
void runge_kutta_4(double a, double b, double y[], int n,
                                  double (*f) (double, double)) {
    double k1, k2, k3, k4, x = a, h = (b - a) / n;

    for (int i = 0; i < n; i++) {
        k1 = f(x, y[i]);
        k2 = f(x + h / 2.0, y[i] + h / 2.0 * k1);
        k3 = f(x + h / 2.0, y[i] + h / 2.0 * k2);
        k4 = f(x + h, y[i] + h * k3);

        y[i+1] = y[i] + h * (k1 + 2.0*k2 + 2.0*k3 + k4) / 6.0;
        x = x + h;
    }
}

// ���ʏo��---------------------------------------------------------
void print_results(double a, double b, double n, double y[]) {
    double h = (b - a) / n;
    int i;

    printf("  x�@        y      \n");
    printf("--------------------\n");
    
    for (i = 0; i <= n; i++) {
        printf("%5.2lf    %8.5lf\n", i*h, y[i]);
    }
}

// ���C���֐�-------------------------------------------------------
int main() {

    double a = 0.0, b = 1.0, y[N+1];

    // �I�C���[�@���s
    y[0] = 1.0;
    euler(a, b, y, N, f);
    printf("�m�I�C���[�@�̌��ʁn\n");
    print_results(a, b, N, y);

    // �z�C���@���s
    y[0] = 1.0;
    heun(a, b, y, N, f);
    printf("\n�m�z�C���@�̌��ʁn\n");
    print_results(a, b, N, y);

    // ���_�@���s
    y[0] = 1.0;
    mid_point(a, b, y, N, f);
    printf("\n�m���_�@�̌��ʁn\n");
    print_results(a, b, N, y);

    // 3�������Q�E�N�b�^�@���s
    y[0] = 1.0;
    runge_kutta_3(a, b, y, N, f);
    printf("\n�m3�������Q�E�N�b�^�@�̌��ʁn\n");
    print_results(a, b, N, y);

    // 4�������Q�E�N�b�^�@���s
    y[0] = 1.0;
    runge_kutta_4(a, b, y, N, f);
    printf("\n�m4�������Q�E�N�b�^�@�̌��ʁn\n");
    print_results(a, b, N, y);

    return 0;
}

