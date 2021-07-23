// Auto_Compensation_Interference.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
#include "Console_String_LIB.h"
#include "Simple_math_LIB.h"
#include "GAUSS_LIB.h"


int main()
{
	double b = 0;//23.45;
	//double b1 = 1.234;
	//double b2 = 5.678;
	//double b3 = 9.10;
	long int ll = 0;
	//char r1[] = "Hello";

	printf("\nb = ");
	scanf("%lf", &b);
	printf("\nll = ");
	scanf("%ld", &ll);
	OU1(b, ll);
	//OU2(ll);
	//OU(b, ll);
	// OF1(b, (int)ll, r1);
	// OFN (b, b1,b2,b3, (int)ll, r1);
	//OFP(b, ll,r1);


	//printf("Hello world!");
	_getch();
    return 0;
}

