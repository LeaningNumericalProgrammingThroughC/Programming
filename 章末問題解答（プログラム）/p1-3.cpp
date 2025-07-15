#include <stdio.h>

int  main() {

	float x, u, err = 0.0, s1 = 0.0, s2 = 0.0;
	int i;

	for (i = 1000000; i >= 1; i--) {

		x = 1.0 / ((float)i * (float)i);
		
		// カハンのアルゴリズム
		u = s1 + (err + x);
		err = x - (u - s1);
		s1 = u;

		// 通常の総和計算
		s2 += x;

	}
	printf("カハンのアルゴリズム ： %.16f\n", s1);
	printf("通常の総和計算　　　 ： %.16f\n", s2);

	return 0;
}