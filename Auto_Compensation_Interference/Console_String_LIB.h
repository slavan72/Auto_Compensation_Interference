# if !defined(__CONSOLE_STRING_LIB_H)
#define __CONSOLE_STRING_LIB_H

/******************************************/
//   Библиотека для составления строк 
//   и вывода этих строк в консоль
//   в определенном порядке
//	 (строим таблицу)
/******************************************/

//#define _USE_MATH_DEFINES
#include "stdafx.h"
//#include <math.h>
//#include <iostream>
//#include <conio.h>
//#include <iostream.h>
//inline double SH(double x){double a;a=exp(x);return (a-1./a)/2.;}
//inline double CH(double x){double a;a=exp(x);return (a+1./a)/2.;}
//*************************************************************************
using namespace std;

/*
Функция OU1(double b, long int ll) производит размещение на экране числа b в формате определяемом
6-ти значным числом ll.
*/
void OU1(double b, long int ll) {//long int ll
								 // b  -выводимое число
								 // ll -6ти-значное число "cdpng" со знаком
	int s, c, p, d; //  n,q
					//float d;
	if (ll<0)s = -1;
	else s = 1;

	ll = s*ll;//если ll<0,то вывод знака (минуса)	 if(ll<0)ll=-1*ll;

			  //парсим  ll -6-ти-значное число "cdpng" 123456. floor() возвращает наибольшее целое число
	c = floor((float)(ll / 10000));          //-количество знаков до запятой = 1234
	d = floor((float)((ll % 10000) / 1000));  //-количество знаков после запятой  = 5
	p = floor((float)((ll % 1000) / 100));     //-количество пробелов после  числа = 6
											   //n=floor((ll%100)/10);        //-количество  getch()
											   //g=ll%10;                     //-количество  переводов строки
	int k, q;//  
	double a = b;
	long int j, l, m, i;  //i

						  //float i;

	if (a<0.) m = -1;
	else m = 1;

	a = a*m;
	l = pow(10, (float)d);//возведение в степень 10 в степени d	(double)(pow(10,d))
	a = a + 0.501 / l;
	i = floor(a); //округление до ближайшего меньшего
	j = floor((a - i)*l);

	if (i == 0) {
		if (d>0)k = c + 1;
		else k = c;
	}
	else {
		q = floor(log((float)i) / log(10.0) + 1);
		k = c - q + 1;
	}

	//сначала выводим пробелы, затем знак, потом целую часть, точку, дробную часть
	for (q = 1; q<k; q++) cout << " ";

	//знак
	if (s<0) {
		if (m<0) cout << "-";
		else cout << " ";
	}

	if (d == 0 || i>0) cout << i;//целая часть   
	if (d>0) cout << ".";//точка
	if (c == 0) goto ZA;

	m = pow(10, (float)d); //
	a = (double)j / m;
	l = 10;

	for (q = 0;q<d;q++) {
		a = a*l;
		k = floor(a + 1. / m);
		cout << k;//дробная часть
		a = a - k;
	}

ZA: for (q = 0; q<p; q++)cout << " ";
} //OU1

  /*
  Функция OU2(long int ll) смещает на экране каретку на порядок определяемый 6-ти значным числом ll.
  */
void OU2(long int ll)
{
	int n, g, q;
	if (ll<0)  ll = -ll; //ll = 123456

	n = floor((float)((ll % 100) / 10));        //n=5  количество  getch() 
	g = ll % 10;                     //g=5 количество  переводов строки 

	for (q = 0; q<n; q++)  _getch();	  // q=6	  { cout<<"a"; getch();} 
	for (q = 0; q<g; q++) cout << "\n"; // q=6	  { cout<<"b"; cout<<"\n";} 
} //OU2

  //***************************************************************************
  /*
  Функции OU размещают на экране double-числа в формате определяемом 6-ти значным числом ll.
  */
void  OU(double a, long int ll) {
	OU1(a, ll);
	OU2(ll);
}

void  OU(double a, double b, long int ll) {
	OU1(a, ll);
	OU1(b, ll);
	OU2(ll);
}
void  OU(double a, double b, double c, long int ll) {
	OU1(a, ll);
	OU1(b, ll);
	OU1(c, ll);
	OU2(ll);
}
void  OU(double a, double b, double c, double d, long int ll) {
	OU1(a, ll);
	OU1(b, ll);
	OU1(c, ll);
	OU1(d, ll);
	OU2(ll);
}
void  OU(double a, double b, double c, double d, double e, long int ll) {
	OU1(a, ll);
	OU1(b, ll);
	OU1(c, ll);
	OU1(d, ll);
	OU1(e, ll);
	OU2(ll);
}
void  OU(double a, double b, double c, double d, double e, double h, long int ll) {
	OU1(a, ll);
	OU1(b, ll);
	OU1(c, ll);
	OU1(d, ll);
	OU1(e, ll);
	OU1(h, ll);
	OU2(ll);
}
//**************************************************************************
/*
ФункциЯ OF1(double b,int ll,char r1[]) производит запись числа b в формате определЯемом
3-x значным числом ll(сo знаком)в строку r1[], НЕ УНИЧТОЖАЯ записанной в нее ранее информации.
*/
void OF1(double b, int ll, char r1[]) {//int ll
									   // b  -выводимое число
									   // ll -3х-значное
	int nc, nd, pr, s;//nd,

	if (ll<0) s = -1;
	else s = 1;

	ll = s*ll;//если ll<0,то вывод знака

	nc = floor((float)(ll / 100));          //количество знаков до запятой
											//if(nc==0)nc=1;
	nd = floor((float)((ll - 100 * nc) / 10));  //количество знаков после запятой
	pr = ll - 100 * nc - 10 * nd;        //количество пробелов после  числа
										 // r1[]-записанное в символьном виде число
										 //char r2[2],r3[2],r4[2],r5[18];
	char r2[4], r3[4], r4[4], r5[48];
	int k, q;
	double a = b;
	long int j, l, m, i; //i,


						 //   strcpy(r1,"");
	strcpy_s(r3, ".");  //strcpy - копировать  "." в r3  
	strcpy_s(r4, " "); //пробел

	if (a<0.) m = -1;
	else m = 1;

	a = a*m;
	l = pow(10, (float)nd); //возведение в степень 10 в степени d  (double)pow(10,nd) 
	a = a + .501 / l;

	//разместить знак в строке
	if (m<0)strcpy(r2, "-");
	else strcpy(r2, " ");

	i = floor(a);	//целая часть
	j = floor((a - i)*l);  //дробная часть

	if (i == 0) {
		if (nd>0)k = nc;
		else k = nc - 1;
	}
	else {
		q = floor(log10((float)i) + 1); k = nc - q;
	}
	for (q = 0; q<k; q++) strcat(r1, r4);

	if (s<0) strcat(r1, r2);//вывод знака
	if (nd == 0 || i >= 0) {
		_itoa(i, r5, 10); //itoa - конвертирует число в строку	 i в десятичном формате в строку r5
		strcat(r1, r5); // к названию добавляем целую часть числа.... strcat - соединяет строки
	}//выводит 0 целых перед точкой
	 //if(nd==0||i>0){itoa(i,r5,10);strcat(r1,r5);}//не выводит 0 целых перед точкой

	if (nd>0) strcat(r1, r3);//затем точку...
	if (nc == 0) goto ZA;

	m = pow(10, (float)nd);  //возведение в степень 10 в степени d (double)pow(10,nd);
	a = (double)j / m;
	l = 10;

	for (q = 0; q<nd; q++) {
		a = a*l;
		k = floor(a + 1. / m);
		_itoa(k, r5, 10);
		strcat(r1, r5);
		a = a - k;
	}
	//printf(r1);//для отладки

ZA:for (q = 0; q<pr; q++)  strcat(r1, r4);
} //OF1

  //***************************************************************************
  /*
  Функции OFN(double a,int ll,char r1[]),..., OFN(double,double,double,double,double,int,char[])
  производЯт запись от одного до пЯти чисел в формате определЯемом 3-x значным числом ll в строку r1[],
  УНИЧТОЖАЯ записанную в нее ранее информацию.
  При работе используют функцию OF1.
  */
void OFN(double a, int ll, char r1[]) {
	strcpy(r1, "");
	OF1(a, ll, r1);

}

void OFN(double a, double b, int ll, char r1[]) {
	strcpy(r1, "");
	OF1(a, ll, r1);
	OF1(b, ll, r1);
}
void OFN(double a, double b, double c, int ll, char r1[]) {
	strcpy(r1, "");
	OF1(a, ll, r1);
	OF1(b, ll, r1);
	OF1(c, ll, r1);
}
void OFN(double a, double b, double c, double d, int ll, char r1[]) {
	strcpy(r1, "");
	OF1(a, ll, r1);
	OF1(b, ll, r1);
	OF1(c, ll, r1);
	OF1(d, ll, r1);
}
void OFN(double a, double b, double c, double d, double e, int ll, char r1[]) {
	strcpy(r1, "");
	OF1(a, ll, r1);
	OF1(b, ll, r1);
	OF1(c, ll, r1);
	OF1(d, ll, r1);
	OF1(e, ll, r1);
}
//***************************************************************************
/*
Функции OFP(double a,int ll,char r1[]),..., OFN(double,double,double,double,double,long int,char[])
производит запись от одного до пЯти чисел в формате определЯемом 3-x значным числом ll в строку r1[],
НЕ УНИЧТОЖАЯ записанной в нее ранее информации.
При работе используют функцию OF1.
*/
void OFP(double a, int ll, char r1[]) {
	OF1(a, ll, r1);
}
void OFP(double a, double b, int ll, char r1[]) {
	OF1(a, ll, r1);
	OF1(b, ll, r1);
}
void OFP(double a, double b, double c, int ll, char r1[]) {
	OF1(a, ll, r1);
	OF1(b, ll, r1);
	OF1(c, ll, r1);
}
void OFP(double a, double b, double c, double d, int ll, char r1[]) {
	OF1(a, ll, r1);
	OF1(b, ll, r1);
	OF1(c, ll, r1);
	OF1(d, ll, r1);
}
void OFP(double a, double b, double c, double d, double e, int ll, char r1[]) {
	OF1(a, ll, r1)
		;OF1(b, ll, r1);
	OF1(c, ll, r1);
	OF1(d, ll, r1);
	OF1(e, ll, r1);
}
//**************************************************************************
#endif __CONSOLE_STRING_LIB_H
