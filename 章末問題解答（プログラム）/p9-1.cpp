#include <stdio.h>
#include <math.h>

#define N 4

void gershgorin_theorem(double a[][N]) {

    int i, j, flag=0;
    double r[N];

    // 半径rを求める
    for (i = 0; i < N; i++) {
        r[i] = 0.0;
        for (j = 0; j < N; j++) {
            if (i != j) r[i] += fabs(a[i][j]);
        }
        printf("R_%d: 中心 %5.2f，  半径 %5.2f\n", i, a[i][i], r[i]);
    }

    // 各円の重なりをチェックする
    printf("\n円の重なりのチェック\n");
    for (i = 0; i < N - 1; i++) {
        for (j = i + 1; j < N; j++) {
            if (fabs(a[i][i] - a[j][j]) <= fabs(r[i] + r[j])) {
                printf("R_%dとR_%dは重なりを持ちます．\n", i, j);
                flag = 1;
            }
        }
    }
    if (flag == 0) {
        printf("各円は重なりを持ちません\n");
    }
}

int main() {

    double a[N][N] = {
        {  2.0, -1.0,  3.0, -4.0 },
        { -2.0,  7.0,  0.0,  3.0 },
        {  0.0,  2.0,  5.0,  3.0 },
        {  2.0,  1.0, -2.0, -6.0 },
    };

    gershgorin_theorem(a);

	return 0;
}