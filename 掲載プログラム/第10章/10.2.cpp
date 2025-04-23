#include <stdio.h>
#include <math.h>

#define N    5           // �s��̎���
#define EPS   1e-8       // ��������p臒l

// N(c)�̌v�Z-------------------------------------------------------------
int n(const double alpha[], double beta[], double c) {

	double w, p_k, p_k_2 = 1.0, p_k_1 = c - alpha[0];
	int i, count = 0;

	if (p_k_1 < 0) count++;
	if (p_k_1 != 0.0) w = p_k_1;

	for (i = 1; i < N; i++) {

		p_k = (c - alpha[i])*p_k_1 - beta[i-1]*beta[i-1]*p_k_2;

		if (p_k * w < 0) count++;
		if (p_k != 0.0) w = p_k;

		p_k_2 = p_k_1;
		p_k_1 = p_k;
	}

	return count;
}

// 2���@------------------------------------------------------------------
double bisection_eigenvalue(int index, double alpha[], double beta[],
	double* min, double* max) {

	double c, lower = *min, upper = *max;

	while ( upper - lower > EPS ) {

		c = (upper + lower) / 2.0;
		if (n(alpha, beta, c) >= index) {
			lower = c;
		} else {
			upper = c;
		}
	}
	
	/* �߂�l���Z�b�g */
	*min = lower; *max = upper;
	return (upper + lower) / 2.0;
}

// ���C���֐�-------------------------------------------------------------
int main() {

	double alpha[N]={ 2,-2,4,1,3 }, beta[N-1]={ 6,-1,9,-5 },c,w;
	double min_0, min, max, lambda[N];
	int i;	

	// ������Ԃ̐ݒ�----------------------------------------------
	min = alpha[0] - fabs(beta[0]);
	max = alpha[0] + fabs(beta[0]);

	for (i = 1; i < N - 1; i++) {
		w = fabs(beta[i - 1]) + fabs(beta[i]);
		if ((alpha[i] - w) < min) {
			min = alpha[i] - w;
		}
		if ((alpha[i] + w) > max) {
			max = alpha[i] + w;
		}
	}
	min_0 = min;

	// �S�Ă̌ŗL�l��2���@�ŋ��߂�----------------------------------
	for (i = 1; i <= N; i++) {

		lambda[i - 1] = bisection_eigenvalue(i, alpha, beta, 
			&min, &max);
		max = min; min = min_0;
	}

	// ���ʏo��-----------------------------------------------------
	printf("�ŗL�l�͈ȉ��̒ʂ�D\n");
	for (i = 0; i < N; i++) {
		printf("%10.6lf\n", lambda[i]);
	}
}