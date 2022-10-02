#include <stdio.h>
#include "../include/math/Zn.hpp"

int main()
{
	Zn<5> num(9);
	int a = 3;
	printf("num = %d\n", num.Value);
	printf("-num = %d\n", (-num).Value);
	printf("num++ = %d\n", (num++).Value);
	printf("num + a = %d\n", (num + a).Value);
	printf("num - a = %d\n", (num - a).Value);
	printf("num * a = %d\n", (num * a).Value);

	return 0;
}
