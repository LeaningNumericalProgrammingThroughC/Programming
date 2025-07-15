#include <stdio.h>

int main() {
	double s = 0.0, x = 0.0, h = 0.7;
	int i;

	for( i = 0; i<=10; i++){

		x = 0.7 * i;
		s += x * x + 1.0;
		x += h;
	}

	printf("s=%lf\n", s);

	return 0;
}