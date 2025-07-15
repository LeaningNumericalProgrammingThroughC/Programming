#include <stdio.h>

int main() {

	double x, a, b, s = 0.0, t = s + 1.0;

	// s+1の演算が正しく出来なくなるまで最大値の探索を続ける
	while (t - s == 1.0) {

		x = 1;
		a = s + x;
		b = a + 1;
		while (b - a == 1.0) {

			// xを倍々に増やしながらs+xが正しく計算できるか確認する
			x = 2.0 * x;
			a = s + x;
			b = a + 1.0;
		}

		// 正しい計算が確認できた増分xを実際にsに足し込む
		s += x / 2.0;
		t = s + 1.0;
	}

	s -= 1.0;
	printf("最大値 s  = %lf\n", s);

	t = s + 1.0;
	printf("t = s + 1 = %lf\n", t);

	printf("t + 1     = %lf\n", t + 1.0);

	return 0;
}