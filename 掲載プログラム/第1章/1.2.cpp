#include <stdio.h>
#include <math.h>

#define N 1000000

int  main() {

	float x, err, u, s1, s2;
	int i;

	err = 0; s1 = 0; s2 = 0;
	for (i = N; i >= 1; i--) {

		x = 1.0 / sqrt((float)i);

		// カハンのアルゴリズム
		u = s1 + (err + x);
		err = x - (u - s1);
		s1 = u;

		// 通常の総和計算
		s2 += x;
	}

	printf("カハンのアルゴリズム ： %.10f\n", s1);
	printf("通常の総和計算　　　 ： %.10f\n", s2);

	return 0;
}