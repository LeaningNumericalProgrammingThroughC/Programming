#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define EPS  1e-8         /* ��������p */
#define MAX_ITER  100     /* �����񐔏�� */

// �֐��̒�`----------------------------------------------------------
double f(double x) {
    return  pow(x, 6.0) + 3 * pow(x, 4.0) - pow(x, 3.0) - 7.0;
}

// 2���@---------------------------------------------------------------
double bisection(double a, double b, int* r) {

    double c, f_a, f_b, f_c;
    int iter = 0;
    
    // �����֐��l�̍Čv�Z������邽�ߕϐ��ɕۑ����Ċ��p����
    f_a = f(a); f_b = f(b); 

    // �����`�F�b�N-------------------------------------------
    if (f_a == 0.0) { *r = iter; return a; }
    if (f_b == 0.0) { *r = iter; return b; }
    if (f_a * f_b > 0.0) {
        printf("���̑��݂𔻒�ł��܂���D\n");
        exit(0); 
    }

    // ��������-----------------------------------------------
    while (fabs(b - a) > EPS) {
        
        iter++;
        c = (a + b) / 2.0;
        f_c = f(c);

        if (f_c == 0.0) break;

        if (f_a * f_c < 0.0) {
            b = c;
        } else {
            a = c; 
            f_a = f_c;
        }
    }

    *r = iter;
    return (a + b) / 2.0;
}

// ���C���֐�----------------------------------------------------------
int main() {

    double a = 0.0, b = 5.0; // �������[a,b]�̐ݒ�
    int iter;

    double c = bisection(a, b, &iter);

    printf("�����񐔂� %d ��C���� %lf �ł��D\n", iter, c);
    
    return 0;
}

