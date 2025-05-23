#include <stdio.h>

int main() {

	double s1, s2;
	int i, flag;

	int N = 1000000000;

	// 第1項から加算
	printf("第1項から加算\n");

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
		printf("  i = %dで加算されなくなりました．　S = %.20lf\n\n", i, s2);
	} else {
		printf("  部分和の計算が完了しました．　S = %.20lf\n\n", s2);
	}

	// 第n項から加算
	printf("第n項から加算\n");

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
		printf("  i = %dで加算されなくなりました．　S = %.20lf\n", i, s2);
	} else {
		printf("  部分和の計算が完了しました．　S = %.20lf\n", s2);
	}

	return 0;
}