# if !defined(__CONSOLE_STRING_LIB_H)
#define __CONSOLE_STRING_LIB_H

/****************************************/
//   ���������� ��� ������ �� �������� 
//   � ������ ���� ����� � �������
//   � ������������ �������
//	 (������ �������)
/****************************************/

//#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <conio.h>
//#include <iostream.h>
//inline double SH(double x){double a;a=exp(x);return (a-1./a)/2.;}
//inline double CH(double x){double a;a=exp(x);return (a+1./a)/2.;}
//*************************************************************************
using namespace std;

void OU1(double b, long int ll) {//long int ll
								 // b  -��������� �����
								 // ll -5��-������� ����� "cdpng" �� ������
	int s, c, p, d; //  n,q
					//float d;
	if (ll<0)s = -1;
	else s = 1;

	ll = s*ll;//���� ll<0,�� ����� ����� (������)

	c = floor((float)(ll / 10000));          //-���������� ������ �� �������. floor() ���������� ���������� ����� �����
	d = floor((float)((ll % 10000) / 1000));  //-���������� ������ ����� �������
	p = floor((float)((ll % 1000) / 100));     //-���������� �������� �����  �����
											   //n=floor((ll%100)/10);        //-����������  getch()
											   //g=ll%10;                     //-����������  ��������� ������
	int k, q;//  
	double a = b;
	long int j, l, m, i;  //i

						  //float i;

	if (a<0.) m = -1;
	else m = 1;

	a = a*m;
	l = pow(10, (float)d);//���������� � ������� 10 � ������� d	(double)(pow(10,d))
	a = a + 0.501 / l;
	i = floor(a); //���������� �� ���������� ��������
	j = floor((a - i)*l);

	if (i == 0) {
		if (d>0)k = c + 1;
		else k = c;
	}
	else {
		q = floor(log((float)i) / log(10.0) + 1);
		k = c - q + 1;
	}

	for (q = 1; q<k; q++) cout << " ";

	//����
	if (s<0) {
		if (m<0) cout << "-";
		else cout << " ";
	}

	if (d == 0 || i>0) cout << i;
	if (d>0) cout << ".";
	if (c == 0) goto ZA;

	m = pow(10, (float)d); //
	a = (double)j / m;
	l = 10;

	for (q = 0;q<d;q++) {
		a = a*l;
		k = floor(a + 1. / m);
		cout << k;
		a = a - k;
	}

ZA: for (q = 0; q<p; q++)cout << " ";
} //OU1

void OU2(long int ll)
{
	int n, g, q;
	if (ll<0)  ll = -ll;

	n = floor((float)((ll % 100) / 10));        //-����������  getch()
	g = ll % 10;                     //-����������  ��������� ������
	for (q = 0; q<n; q++) _getch();
	for (q = 0; q<g; q++) cout << "\n";
} //OU2

  //***************************************************************************
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
������� OF1(double b,int ll,char r1[]) ���������� ������ ����� b � ������� ������������
3-x ������� ������ ll(�o ������)� ������ r1[], �� ��������� ���������� � ��� ����� ����������.
*/
void OF1(double b, int ll, char r1[]) {//int ll
									   // b  -��������� �����
									   // ll -3�-�������
	int nc, nd, pr, s;//nd,

	if (ll<0) s = -1;
	else s = 1;

	ll = s*ll;//���� ll<0,�� ����� �����

	nc = floor((float)(ll / 100));          //-���������� ������ �� �������
											//if(nc==0)nc=1;
	nd = floor((float)((ll - 100 * nc) / 10));  //-���������� ������ ����� �������
	pr = ll - 100 * nc - 10 * nd;        //-���������� �������� �����  �����
										 // r1[]-���������� � ���������� ���� �����
										 //char r2[2],r3[2],r4[2],r5[18];
	char r2[4], r3[4], r4[4], r5[48];
	int k, q;
	double a = b;
	long int j, l, m, i; //i,

						 //   strcpy(r1,"");
	strcpy(r3, ".");
	strcpy(r4, " "); //������

	if (a<0.) m = -1;
	else m = 1;

	a = a*m;
	l = pow(10, (float)nd); //���������� � ������� 10 � ������� d  (double)pow(10,nd) 
	a = a + .501 / l;

	if (m<0)strcpy(r2, "-");
	else strcpy(r2, " "); //����

	i = floor(a);
	j = floor((a - i)*l);

	if (i == 0) {
		if (nd>0)k = nc;
		else k = nc - 1;
	}
	else {
		q = floor(log10((float)i) + 1); k = nc - q;
	}
	for (q = 0; q<k; q++) strcat(r1, r4);

	if (s<0) strcat(r1, r2);//����� �����
	if (nd == 0 || i >= 0) {
		itoa(i, r5, 10);strcat(r1, r5);
	}//������� 0 ����� ����� ������
	 //if(nd==0||i>0){itoa(i,r5,10);strcat(r1,r5);}//�� ������� 0 ����� ����� ������

	if (nd>0) strcat(r1, r3);
	if (nc == 0) goto ZA;

	m = pow(10, (float)nd);  //���������� � ������� 10 � ������� d (double)pow(10,nd);
	a = (double)j / m;
	l = 10;

	for (q = 0; q<nd; q++) {
		a = a*l;
		k = floor(a + 1. / m);
		itoa(k, r5, 10);
		strcat(r1, r5);
		a = a - k;
	}
ZA:for (q = 0; q<pr; q++)  strcat(r1, r4);
} //OF1

  //***************************************************************************
  /*
  ������� OFN(double a,int ll,char r1[]),...,
  OFN(double,double,double,double,double,int,char[]) ����������
  ������ �� ������ �� ���� ����� � ������� ������������ 3-x �������
  ������ ll � ������ r1[],
  ���������
  ���������� � ��� ����� ����������.
  ��� ������ ���������� ������� OF1.
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
������� OFP(double a,int ll,char r1[]),...,
OFN(double,double,double,double,double,long int,char[]) ����������
������ �� ������ �� ���� ����� � ������� ������������ 3-x �������
������ ll � ������ r1[],
�� ���������
���������� � ��� ����� ����������.
��� ������ ���������� ������� OF1.
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
