#include <stdio.h>

int main() {

	int a[4] = { 2,4,3,-1 }, b[4] = { 5,-2,2,1 };
	int i, j, k, p, s, f;

	p = 1;
	// Piのループ
	for (i = 0; i <= 3; i++) {
		
		s = 0;
		// Sumのループ
		for (j = 0; j <= 3; j++) {
		
			f = 1;
			// (j+1)!のループ
			for (k = 1; k <= j + 1; k++) {
				f *= k;
			}

			s += f * a[i] * b[j];
		}

		f = 1;
		// (i+2)!のループ
		for (j = 1; j <= i+2; j++) {
			f *= j;
		}

		p *= s + f;
	}

	printf("計算結果：%d\n", p);

	return 0;
}