#include <stdio.h>

int main() {

	double s1, s2;
	int i, flag;

	int N = 1000000000;

	// ‘æ1€‚©‚ç‰ÁZ
	printf("‘æ1€‚©‚ç‰ÁZ\n");

	flag = 0;
	s1 = s2 = 0;
	for (i = 1; i <= N; i++) {

		s2 = s1 + 1.0 / ((double)i * (double)i);

		if (s1 == s2) {
			flag = 1;
			break;
		}

		s1 = s2;
	}

	if (flag == 1) {
		printf("  i = %d‚Å‰ÁZ‚³‚ê‚È‚­‚È‚è‚Ü‚µ‚½D@S = %.20lf\n\n", i, s2);
	} else {
		printf("  •”•ª˜a‚ÌŒvZ‚ªŠ®—¹‚µ‚Ü‚µ‚½D@S = %.20lf\n\n", s2);
	}

	// ‘æn€‚©‚ç‰ÁZ
	printf("‘æn€‚©‚ç‰ÁZ\n");

	flag = 0;
	s1 = s2 = 0;
	for (i = N; i >= 1; i--) {

		s2 = s1 + 1.0 / ((double)i * (double)i);

		if (s1 == s2) {
			flag = 1;
			break;
		}

		s1 = s2;
	}

	if (flag == 1) {
		printf("  i = %d‚Å‰ÁZ‚³‚ê‚È‚­‚È‚è‚Ü‚µ‚½D@S = %.20lf\n", i, s2);
	} else {
		printf("  •”•ª˜a‚ÌŒvZ‚ªŠ®—¹‚µ‚Ü‚µ‚½D@S = %.20lf\n", s2);
	}

	return 0;
}