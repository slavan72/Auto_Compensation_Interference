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

//********************************************************
void UPOMCAN(double dt, double fo, long int n, int N, double zad, double v[2], double af[2], double pole[2], double&u) {//         v-  2 аргумента для RALLEY
 //        af- случайная амп и фаза oгибающей
 //выход:	u-  мгновенное напряжение
	int m;

	m = (int)floor(zad / dt);
	if (n%N == m + 1) {
		RALLEY(v, af);
		if (n < N)  RALLEY(v, af);
	}

	u = af[0] * pole[0] * sin(2 * M_PI*fo*(n*dt/*-zad*/) + pole[1] + af[1]);

}//UPOMCAN

 //********************************************************
void SHUM(double dt, double fo, long int n, int Nsl, double tns,
	double vs[2], double afs[2], double&u)
{ //u- нормир шум
	int m;

	m = (int)floor(tns*Nsl);
	//cout<<"m+1="<<m+1<<"  n%Nsl="<<n%Nsl<<"\n";
	if (n%Nsl == (m + 1) % Nsl) {
		RALLEY(vs, afs); if (n < Nsl)RALLEY(vs, afs);
		//cout<<"\ntjrop\n";
	}

	u = afs[0] * sin(2 * M_PI*fo*(n - m)*dt + afs[1]);

}//SHUM

 //********************************************************
void SIGNAL(double dt, double fo, long int n, int Nsl, double kalibr, int Nc, int Np, double v[4][5][2], double af[4][5][2],
	double tns[5], double vs[5][2], double afs[5][2], double zad[4][5], double pole[4][5][2], double usignal[5]) {
	int nc, np, m, q;       //double uv_f,u_f;
	int nomkontip = 0;
	//cout<<"shum  ";
	for (nc = 0;nc < Nc;nc++) {
		//cout<<afs[nc][0]<<"  "<<afs[nc][1]<<"  "<<usignal[nc]<<"  ";
		//cout<<"signal_0"<<"\n";
		SHUM(dt, fo, n, Nsl, tns[nc], vs[nc], afs[nc], u_f);
		usignal[nc] = u_f;
		//cout<<nc<<vs[nc][0]<<"  "<<vs[nc][1]<<"  "<<afs[nc][0]<<"  "<<afs[nc][1]<<"  "<<usignal[nc]<<"  ";
		for (np = 0;np < Np;np++) {
			UPOMCAN(dt, fo, n, Nsl, zad[np][nc], v[np][nc], af[np][nc], pole[np][nc], u_f);
			usignal[nc] = usignal[nc] + u_f;
			//cout<<"nc="<<nc<<"  np="<<np<<"  ucan="<<ucan[nc]<<"  zad="<<zad[np][nc]<<"\n";getch();
		}//np
		usignal[nc] = usignal[nc] * ashum*kalibr;

		// UFV(nc,0,FVST,FVCT,FVSN,FVCN,usignal[nc],usignals[nc],usignalc[nc]);
		//  usignals[nc]=2*usignals[nc];
		//  usignalc[nc]=2*usignalc[nc];

		uv_f = usignal[nc];
		//-------------------------------------------------------------
		VHODFIL_2(0, nc, KT, KN, uv_f, u_f);
		usignal[nc] = u_f;
		//-------------------------------------------------------------
		//cout<<usignal[nc]<<" ";
	}//nc
	 //cout<<"\n"; getch();
}//SIGNAL
