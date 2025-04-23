#include <stdio.h>
#include <math.h>

#define N 100            // ������
#define MAX_ITER 10      // �ő唽����
#define EPS 1e-6         // ��������p臒l

// �����������̉E�ӂ��`----------------------------------------------------------
double f(double x, double y) {

	return 8*x*y;
}

// �����Q�E�N�b�^�@�i1�X�e�b�v���s�j-----------------------------------------------
double runge_kutta_1_step(double x, double y, double h, 
	                                 double (*f)(double, double)) {
	double k1, k2, k3, k4;
	
	k1 = f(x, y);
	k2 = f(x + h / 2.0, y + h / 2.0 * k1);
	k3 = f(x + h / 2.0, y + h / 2.0 * k2);
	k4 = f(x + h, y + h * k3);

	return y +  h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

// �\���q�C���q�@------------------------------------------------------------------
void predictor_corrector_method(double a, double b, int n, double y[], 
                               int num_iter[], double (*f)(double, double)) {

	double y_new, y_old, h = (b-a)/n;
	double f_i[5] ;
	int i, j, iter=0;

	// �����Q�E�N�b�^�@�ɂ��ŏ���3�_�����߂�-------------------------------
	f_i[0] = f(a , y[0]);
	for (i = 1; i < 4; i++) {
		y[i] = runge_kutta_1_step(a+(i-1)*h, y[i-1], h, f);
		f_i[i] = f(a + i * h, y[i]);
	}

	// �A�_���X�̕��@��4�_�ڈȍ~�����߂�-------------------------------------
	for (i = 4; i <= n; i++) {

		// �A�_���X�E�o�b�V���t�H�[�X�@---------------------------------
		y_old = y[i-1] + h * (-9.0*f_i[0] + 37.0*f_i[1] 
			                 - 59.0*f_i[2] + 55.0*f_i[3]) / 24.0;
		
		// �A�_���X�E�����g���@�Ŕ�������-------------------------------
		iter = 0;
		do {
			f_i[4] = f(a+i*h, y_old);
			y_new = y[i-1] + h / 24.0 * (f_i[1] - 5.0 * f_i[2] 
				                  + 19.0 * f_i[3] + 9.0 * f_i[4]);
			if (fabs(y_new - y_old) < EPS) break;
			y_old = y_new;
			iter++;
		} while (iter < MAX_ITER);

		// �������Ȃ��ꍇ�̓��b�Z�[�W��\��-----------------------------
		if (iter == MAX_ITER) {
			printf("��%d�X�e�b�v�F�������܂���ł����D\n", i);
		}

		for (j = 0; j < 3; j++)	f_i[j] = f_i[j + 1];
		y[i] = y_new;
		f_i[3] = f(a + i * h, y[i]);

		num_iter[i] = iter;
	}
}

// ���C���֐�----------------------------------------------------------------------
int main() {

	double a = 0, b = 1, y[N+1];
	int i, num_iter[N + 1] = { 0 };

	y[0] = 1.0;
	predictor_corrector_method(a, b, N, y, num_iter, f);

	// ���ʏo��------------------------------------------------------------
	printf("   x         y      ������\n");
	printf("------------------------------\n");
	for (i = 0; i <= N; i++) {
		printf("%5.2lf    %8.5lf      %d\n", (b-a)/N*i, y[i], num_iter[i]);
	}

	return 0;
}