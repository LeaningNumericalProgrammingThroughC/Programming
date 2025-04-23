#include <stdio.h>

int main() {

	double s1, s2, p1, p2;
	int i, j, k;

	s1 = 0.0;
	for (i = 1; i <= 6; i++) {

		p1 = 1.0;
		for (j = 1; j <= i; j++) {
			p1 += 3.0 * j + 1.0;
		}

		s2 = 0.0;
		for (j = 1; j <= i; j++) {

			p2 = 1.0;
			for (k = 1; k < j; k++) {
				p2 *= k;
			}

			s2 += p2 + 3.0;
		}

		s1 += p1 + 2.0 * i * s2;
	}

	printf("ŒvŽZŒ‹‰Ê‚Í%lf‚Å‚·D\n", s1);

	return 0;
}