# if !defined(__GAUSS_LIB_H)
#define __GAUSS_LIB_H

#define _USE_MATH_DEFINES

/******************************************/
// ���������� ������� ��� ������� ������� 
// �������� ��������� � ������������      
// �������������� (GAUSS), � ����� �������  
// ���������� ������� � �o���������� 
// �������������� � ������������ ���� (TREGAUSSC)
//  WLIB3.h
/******************************************/

#include "stdafx.h"

//#include <math.h>
#include <complex>
//using std::complex;

using namespace std;

//****************************************************************
void GAUSS2c(int const I, complex <double> a[], complex <double> x[])
{//GAUSS: ����������� ������ ������� �������� ��������� � ������������
 //       ��������������
 //���������: GAUSS1c(I,A,X)
 //     const  I     - ����������� �������
 //     const  M
 //     double A[]   -����������� ������� �������������
 //            �-      ����� �������
 //cout<<"I="<<I<<"  M="<<M<<"\n";                     getch();
	int i, j, m, L; //,i1,j1
	complex <double> b;
	L = I + 1;
	//---------------------------------------------------------------
	/*
	cout<<" ����������� �������  I="<<I<<"\n\n";
	cout<<"����������� ������� ������������� \n";
	for(m=0;m<I*I+I;m++)
	{cout<<a[m];
	if((m+1)%(I+1)==0)cout<<"\n";
	}
	cout<<"\n\n";
	*/
	//---------------------------------------------------------------

	for (m = 1; m<I; m++) {
		if (a[(m - 1)*I + m - 2 + m] == complex <double>(0.0, 0.0)) {

			for (i = m; i<I; i++) { //i-����� ��-�

				if (a[i*L + m - 1] != complex <double>(0.0, 0.0)) {
				}
			}
		}
		for (i = m; i<I; i++) {

			b = a[(m - 1)*I + m - 2 + m] / a[i*L + m - 1];
			//cout<<a[(m-1)*I+m-2+m]<<"  "<<a[i*L+m-1]<<"  b="<<b<<"\n";
			for (j = m - 1;j<L;j++) {

				a[i*L + j] = a[i*L + j] * b - a[(m - 1)*L + j];

				//for(i1=0;i1<I;i1++)                                       //�����
				//{for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1];cout<<"\n";} //�������������
				//cout<<"\n";getch();                                       //������
			}
		}
	}
	// cout<<"������� ����������� �������\n"; getch();
	for (j = I - 1; j>0; j--) {

		for (i = j - 1; i >= 0; i--) {

			b = a[j*L + I] / a[j*L + j];
			a[i*L + I] = a[i*L + I] - b*a[i*L + j];
			//      cout<<"j="<<j<<"  i="<<i<<"\n";                     //�����
			//      for(i1=0;i1<I;i1++)                                  // |
			//      { for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1];cout<<"\n";}  //�������������
			//      cout<<"\n";getch();                                  //������
		}
	}
	cout << "�����:";

	for (i = 0; i<I; i++) {
		x[i] = a[i*L + I] / a[i*L + i];
		//    cout<< x[i]<<"\n"; getch();
	}//i
	 //cout<<"\n";    getch();
}  //����� GAUSS2c


   //*************************************************************************
   //************************************************************************
void GAUSS1c(int const I, complex <double> a[], complex <double> x[])
{//GAUSS: ����������� ������ ������� �������� ��������� � ������������
 //       ��������������
 //���������: GAUSS1c(I,A,X)
 //     const  I     - ����������� �������
 //     const  M
 //     double A[]   -����������� ������� �������������
 //            �-      ����� �������
 //cout<<"I="<<I<<"  M="<<M<<"\n";                     getch();
	int i, j, m, i1, j1, L;

	complex <double> b;
	L = I + 1;
	for (m = 1; m<I; m++) {
		for (i = m; i<I; i++) {
			if (a[i*L + m - 1] == complex <double>(0.0, 0.0)) goto ZZZ;
			b = a[(m - 1)*I + m - 2 + m] / a[i*L + m - 1];
			//cout<<a[(m-1)*I+m-2+m]<<"  "<<a[i*L+m-1]<<"  b="<<b<<"\n";
			for (j = m - 1;j<L;j++) {
				a[i*L + j] = a[i*L + j] * b - a[(m - 1)*L + j];
				// for(i1=0;i1<I;i1++)                                       //�����
				//  {for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1]<<"  ";cout<<"\n";} //�������������
				//   cout<<"\n";getch();                                       //������
			}
		ZZZ:;
		}
	}  // cout<<"������� ����������� �������\n"; getch();

	for (j = I - 1; j>0; j--) {
		for (i = j - 1; i >= 0; i--) {
			b = a[j*L + I] / a[j*L + j];
			a[i*L + I] = a[i*L + I] - b*a[i*L + j];
			//    cout<<"i="<<i<<"  j="<<j<<"\n";    for(i1=0;i1<I;i1++)     //�����
			//  { for(j1=0;j1<=I;j1++)cout<<a[i1*(I+1)+j1]<<"  ";cout<<"\n";}  //�������������
			//  cout<<"\n";getch();                                       //������
		}
	}
	//cout<<"�����:";
	for (i = 0; i<I; i++) {
		x[i] = a[i*L + I] / a[i*L + i];
		//cout<< x[i]<<"\n"; getch();
	}//i
	 //cout<<"\n";    getch();
}  //����� GAUSS1c


   //*************************************************************************
void GAUSSc(int const I, long double ar[], long double ai[],
	long double xr[], long double xi[], long double det[])
{//GAUSSc: o���������� �o���� �o�������o� ������� �������� ���������
 //���������: GAUSSc(I,Ar,Ai,Xr,Xi)
 //     const  I          - �������o��� �������
 //long double �r[I*I+I]  -����������� ������� ������������o� �����
 //                        �o���������o�
 //long double �i[I*I+I]  -����������� ������� ����o� ����� �o���������o�
 //long double �r[I]      -����o� ������������o� ����� �o���� �������
 //long double Xi[I]      -����o� ����o� ����� �o���� �������
 //long double det[2]     -��������������[0] � ������[1] ����� ������������
 //cout<<"I="<<I<<"\n";                     getch();
	int i, j, m, p, q, i1, j1, L;
	long double br, bi, cr, ci, dr, di, e;
	L = I + 1;
	for (m = 1; m<I; m++) {
		for (i = m; i<I; i++) {
			p = (m - 1)*I + m + m - 2;
			cr = ar[p];
			ci = ai[p];
			q = i*L + m - 1;
			dr = ar[q];
			di = ai[q];
			e = 1 / (dr*dr + di*di);
			br = (cr*dr + ci*di)*e;
			bi = (ci*dr - cr*di)*e;
			// b=a[(m-1)*I+m-2+m]/a[i*L+m-1];
			for (j = m - 1;j<L;j++) {
				p = i*L + j;
				cr = ar[p];
				ci = ai[p];
				q = (m - 1)*L + j;
				ar[p] = cr*br - ci*bi - ar[q];
				ai[p] = ci*br + cr*bi - ai[q];
				//a[i*L+j]=a[i*L+j]*b-a[(m-1)*L+j];
				// for(i1=0;i1<I;i1++)                                       //���o�
				//  {for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1]<<"  ";cout<<"\n";} //��o�����o����
				//  cout<<"\n";getch();                                       //������
			}
		}
	}  // cout<<"������� �����o����� �������\n"; getch();
	   //-----------------------------------------------------------------
	   //        ���������� �������� ������������
	det[0] = det[1] = 1;
	for (i1 = 0; i1<I; i1++)
	{
		p = i1*L + i1; //����� ��������
		cr = ar[p] * det[0] - ai[p] * det[1]; // real
		det[1] = ar[p] * det[1] + ai[p] * det[0];//image
		det[0] = cr;
	}
	//-----------------------------------------------------------------
	for (j = I - 1; j>0; j--) {
		for (i = j - 1; i >= 0; i--) {
			p = j*L + I;
			cr = ar[p];
			ci = ai[p];
			q = j*L + j;
			dr = ar[q];
			di = ai[q];
			e = 1 / (dr*dr + di*di);
			br = (cr*dr + ci*di)*e;
			bi = (ci*dr - cr*di)*e;
			//b=a[j*L+I]/a[j*L+j];
			p = i*L + I; cr = ar[p];
			ci = ai[p];
			q = i*L + j;
			dr = ar[q];
			di = ai[q];
			ar[p] = ar[p] - br*dr + di*bi;
			ai[p] = ai[p] - bi*dr - br*di;
			// a[i*L+I]=a[i*L+I]-b*a[i*L+j];
			//    cout<<"i="<<i<<"  j="<<j<<"\n";    for(i1=0;i1<I;i1++)     //���o�
			//  { for(j1=0;j1<=I;j1++)cout<<a[i1*(I+1)+j1]<<"  ";cout<<"\n";}  //��������o����
			//  cout<<"\n";getch();                                       //������
		}
	}
	//cout<<"�o���:";
	for (i = 0; i<I; i++) {
		p = i*L + I;
		cr = ar[p];
		ci = ai[p];
		q = i*L + i;
		dr = ar[q];
		di = ai[q];
		e = 1 / (dr*dr + di*di);
		xr[i] = (cr*dr + ci*di)*e;
		xi[i] = (ci*dr - cr*di)*e;
		//  cout<< xr[i]<<"  "<<xi[i]<<"\n"; getch();
	}
	//    cout<<"\n";    getch();
}  //�o��� GAUSSc
   //*************************************************************************
void GAUSScc(int I, long double ar[], long double ai[],
	double xr[], double xi[], long double det[])
{//GAUSSc: o���������� �o���� �o�������o� ������� �������� ���������
 //���������: GAUSSc(I,Ar,Ai,Xr,Xi)
 //     const  I          - �������o��� �������
 //long double �r[I*I+I]  -����������� ������� ������������o� �����
 //                        �o���������o�
 //long double �i[I*I+I]  -����������� ������� ����o� ����� �o���������o�
 //long double �r[I]      -����o� ������������o� ����� �o���� �������
 //long double Xi[I]      -����o� ����o� ����� �o���� �������
 //long double det[2]     -��������������[0] � ������[1] ����� ������������
 //cout<<"I="<<I<<"\n";                     getch();
	int i, j, m, p, q, i1, j1, L;
	long double br, bi, cr, ci, dr, di, e;
	L = I + 1;
	for (m = 1; m<I; m++) {
		for (i = m; i<I; i++) {
			p = (m - 1)*I + m + m - 2;
			cr = ar[p];
			ci = ai[p];
			q = i*L + m - 1;
			dr = ar[q];
			di = ai[q];
			e = 1 / (dr*dr + di*di);
			br = (cr*dr + ci*di)*e;
			bi = (ci*dr - cr*di)*e;
			// b=a[(m-1)*I+m-2+m]/a[i*L+m-1];
			for (j = m - 1;j<L;j++) {
				p = i*L + j;
				cr = ar[p];
				ci = ai[p];
				q = (m - 1)*L + j;
				ar[p] = cr*br - ci*bi - ar[q];
				ai[p] = ci*br + cr*bi - ai[q];
				//a[i*L+j]=a[i*L+j]*b-a[(m-1)*L+j];
			}
		}
	} //  cout<<"������� �����o����� �������\n"; getch();

	for (i1 = 0; i1<I; i1++) {
		for (j1 = 0; j1 <= I; j1++)
			cout << ar[i1*L + j1] << "+i" << ai[i1*L + j1] << "  ";   //��o�����o����	 ������
		cout << "\n";
	}
	cout << "\n";
	getch();
	//-----------------------------------------------------------------
	/*
	//        ���������� �������� ������������
	det[0]=det[1]=1;
	for(i1=0;i1<I;i1++)
	{p=i1*L+i1; //����� ��������
	cr=ar[p]*det[0]-ai[p]*det[1]; // real
	det[1]=ar[p]*det[1]+ai[p]*det[0];//image
	det[0]=cr;
	}
	*/
	//-----------------------------------------------------------------
	//        ���������� �������� ������������
	det[0] = ar[0];
	det[1] = ai[0];
	for (i1 = 1; i1<I; i1++) {
		p = i1*L + i1; //����� ��������
		cr = ar[p] * det[0] - ai[p] * det[1]; // real
		det[1] = ar[p] * det[1] + ai[p] * det[0];//image
		det[0] = cr;
	}
	//-----------------------------------------------------------------
	for (j = I - 1; j>0; j--) {
		for (i = j - 1; i >= 0; i--) {
			p = j*L + I;
			cr = ar[p];
			ci = ai[p];
			q = j*L + j;
			dr = ar[q];
			di = ai[q];
			e = 1 / (dr*dr + di*di);
			br = (cr*dr + ci*di)*e;
			bi = (ci*dr - cr*di)*e;
			//b=a[j*L+I]/a[j*L+j];
			p = i*L + I;
			cr = ar[p];
			ci = ai[p];
			q = i*L + j;
			dr = ar[q];
			di = ai[q];
			ar[p] = ar[p] - br*dr + di*bi;
			ai[p] = ai[p] - bi*dr - br*di;
			// a[i*L+I]=a[i*L+I]-b*a[i*L+j];
			//    cout<<"i="<<i<<"  j="<<j<<"\n";    for(i1=0;i1<I;i1++)     //���o�
			//  { for(j1=0;j1<=I;j1++)cout<<a[i1*(I+1)+j1]<<"  ";cout<<"\n";}  //��������o����
			//  cout<<"\n";getch();                                       //������
		}
	}
	//cout<<"�o���:";
	for (i = 0; i<I; i++) {
		p = i*L + I;
		cr = ar[p];
		ci = ai[p];
		q = i*L + i;
		dr = ar[q];
		di = ai[q];
		e = 1 / (dr*dr + di*di);
		xr[i] = (cr*dr + ci*di)*e;
		xi[i] = (ci*dr - cr*di)*e;
		//  cout<< xr[i]<<"  "<<xi[i]<<"\n"; getch();
	}
	//    cout<<"\n";    getch();
}  //�o��� GAUSScc


   //*************************************************************************
   //*************************************************************************
void GAUSSccc(int I, long double ar[], long double ai[],
	double xr[], double xi[], long double det[])
{//GAUSSc: o���������� �o���� �o�������o� ������� �������� ���������
 //���������: GAUSSc(I,Ar,Ai,Xr,Xi)
 //             ������� ar,ai �� �����������
 //     const  I          - �������o��� �������
 //long double �r[I*I+I]  -����������� ������� ������������o� �����
 //                        �o���������o�
 //long double �i[I*I+I]  -����������� ������� ����o� ����� �o���������o�
 //long double �r[I]      -����o� ������������o� ����� �o���� �������
 //long double Xi[I]      -����o� ����o� ����� �o���� �������
 //long double det[2]     -��������������[0] � ������[1] ����� ������������
 //cout<<"I="<<I<<"\n";                     getch();

	double fr[30], fi[30], ftr[30], fti[30];
	// ofstream DTmat("MATRIX.dat");
	int i, j, m, p, q, s, t, i1, j1, L;
	long double br, bi, cr, ci, dr, di, e;
	L = I + 1;
	if (I*L>30) {
		cout << " ����������\n"
			"��������� \n"
			" ����� �������� fr,fi\n"
			" � ��������� GAUSSccc\n";
		getch();
		getch();
	}

	//DTmat<<"�������� �������\n";
	for (i1 = 0; i1<I; i1++) {                                       //���o�
		for (j1 = 0; j1 <= I; j1++) {
			//cout<<i1<<"  "<<j1<<"  "<<i1*L+j1<<"\n"; getch();
			fr[i1*L + j1] = ar[i1*L + j1];fi[i1*L + j1] = ai[i1*L + j1];
			//DTmat<<ar[i1*L+j1]<<"+i"<<ai[i1*L+j1]<<"  ";
		}
		//    DTmat<<"\n";
	}                                           //��o�����o����
												//  DTmat<<"\n";
	for (m = 0; m<I; m++) {   //����� ����������� �������
		q = (I + 2)*m;//����� ������������ ��-�� ����������� �������
		e = 1 / (ar[q] * ar[q] + ai[q] * ai[q]);
		for (i = m + 1;i<I;i++) { //����� ������������� ������
			s = i*L + m;           //����� ��-�� [i,m];
			br = (ar[q] * ar[s] + ai[q] * ai[s])*e;
			bi = -(ai[q] * ar[s] - ar[q] * ai[s])*e;
			//��������� ��� ���������� ������
			for (j = m; j<L; j++) {  //����� �������������� �������
				p = i*L + j;         //����� �������������� ��������
				t = p - L*(i - m);               //����� ����������� ��������
				ar[p] = ar[p] - (ar[t] * br - ai[t] * bi);
				if (fabs(ar[p])<1e-17)ar[p] = 0;
				ai[p] = ai[p] - (ai[t] * br + ar[t] * bi);
				if (fabs(ai[p])<1e-17)ai[p] = 0;
				//cout<<"m="<<m<<"  i="<<i<<"  j="<<j<<"  q="<<q<<"   s="<<s<<"   p="<<p<<
				//      "  ar="<<ar[p]<<"  ai="<<ai[p]<< "\n";
			}
		}
	} //  cout<<"������� �����o����� �������\n";
	  //getch();
	  //DTmat<<"����������� �������\n";
	for (i1 = 0; i1<I; i1++) {                                       //���o�
		for (j1 = 0; j1 <= I; j1++) {
			//DTmat<<ar[i1*L+j1]<<"+i"<<ai[i1*L+j1]<<"  ";
			ftr[i1*L + j1] = ar[i1*L + j1];
			fti[i1*L + j1] = ai[i1*L + j1];
		}
		//DTmat<<"\n";
	}                                           //��o�����o����
												//DTmat<<"\n";
												//getch();                                       //������
												//-----------------------------------------------------------------
												/*
												//        ���������� �������� ������������
												det[0]=det[1]=1;
												for(i1=0;i1<I;i1++)
												{p=i1*L+i1; //����� ��������
												cr=ar[p]*det[0]-ai[p]*det[1]; // real
												det[1]=ar[p]*det[1]+ai[p]*det[0];//image
												det[0]=cr;
												}
												*/
												//-----------------------------------------------------------------
												//        ���������� �������� ������������
	det[0] = ar[0];det[1] = ai[0];
	for (i1 = 1; i1<I; i1++) {
		p = i1*L + i1; //����� ��������
		cr = ar[p] * det[0] - ai[p] * det[1]; // real
		det[1] = ar[p] * det[1] + ai[p] * det[0];//image
		det[0] = cr;
	}
	//-----------------------------------------------------------------
	for (m = I - 1; m>0; m--) {      //����� ����������� �������
		q = (I + 2)*m;//����� ������������ ��-�� ����������� �������
		e = 1 / (ar[q] * ar[q] + ai[q] * ai[q]);

		for (i = m - 1; i >= 0; i--) { //����� ������������� ������
			s = i*L + m;           //����� ��-�� [i,m];
			br = (ar[q] * ar[s] + ai[q] * ai[s])*e;
			bi = -(ai[q] * ar[s] - ar[q] * ai[s])*e;
			//��������� ��� ���������� ������
			p = q - (m - i)*L;     //����� �������������� �������� q-�� �������

			ar[p] = ar[p] - (ar[q] * br - ai[q] * bi);
			if (fabs(ar[p])<1e-17) ar[p] = 0;

			ai[p] = ai[p] - (ai[q] * br + ar[q] * bi);
			if (fabs(ai[p])<1e-17) ai[p] = 0;

			t = L*(i + 1) - 1;   //����� �������������� �������� I-�� �������
			t = p + I - m;

			ar[t] = ar[t] - (ar[t + L*(m - i)] * br - ai[t + L*(m - i)] * bi);
			if (fabs(ar[t])<1e-17)  ar[t] = 0;

			ai[t] = ai[t] - (ai[t + L*(m - i)] * br + ar[t + L*(m - i)] * bi);
			if (fabs(ai[t])<1e-17)  ai[t] = 0;
			//------------------------------------------------------------------------
			//            ������������� ������
			/*
			DTmat<<"�������  m="<<m<<"\n���� �� q="<<q<<" "<<ar[q]<<"+i"<<ai[q]
			<<"\n������ i="<<i<<"\n��� s="<<s<<"  "<<ar[s]<<"+i"<<ai[s]
			<<"\n������ ��� p="<<p<<"  "<<br<<"+i"<<bi
			<<"\n���������� ��� q="<<q<<"  "<<ar[q]<<"+i"<<ai[q]
			<<"\n������ ��� p="<<p<<"  "<<ar[p]<<"+i"<<ai[p]
			<<"\n���������� ��� t+L*(m-i)="<<t+L*(m-i)<<"  "
			<<ar[t+L*(m-i)]<<"+i"<<ai[t+L*(m-i)];
			DTmat<<"\n������ ��� t="<<t<<"  "<<ar[t]<<"+i"<<ai[t]<<"\n\n";
			// a[i*L+I]=a[i*L+I]-b*a[i*L+j];
			//    cout<<"i="<<i<<"  j="<<j<<"\n";    for(i1=0;i1<I;i1++)     //���o�
			//  { for(j1=0;j1<=I;j1++)cout<<a[i1*(I+1)+j1]<<"  ";cout<<"\n";}  //��������o����
			//  cout<<"\n";getch();                                       //������
			DTmat<<"������������� ���������� ������� \n";
			for(i1=0;i1<I;i1++)                                       //���o�
			{for(j1=0;j1<=I;j1++)DTmat<<ar[i1*L+j1]<<"+i"<<ai[i1*L+j1]<<"  ";
			DTmat<<"\n";
			}                                           //��o�����o����
			DTmat<<"\n";
			//                         */
			//---------------------------------------------------------------------------
		}
	}
	//cout<<"�o���:";
	for (i = 0; i<I; i++) {
		p = i*L + I;
		cr = ar[p];
		ci = ai[p];
		q = i*L + i;
		dr = ar[q];
		di = ai[q];
		e = 1 / (dr*dr + di*di);
		xr[i] = (cr*dr + ci*di)*e;
		xi[i] = (ci*dr - cr*di)*e;
		//  cout<< xr[i]<<"  "<<xi[i]<<"\n"; getch();
	}
	//DTmat<<"���������� �������                                         �����\n";
	for (i1 = 0; i1<I; i1++) {                                       //���o�
		for (j1 = 0; j1 <= I; j1++);//DTmat<<ar[i1*L+j1]<<"+i"<<ai[i1*L+j1]<<"  ";
									//DTmat<<"         "<<xr[i1]<<"+i"<<xi[i1]<<"\n"

	}                                           //��o�����o����
												//  DTmat<<"\n";
												//DTmat<<"�������� � �������� �������\n";
	for (i1 = 0; i1<I; i1++) {        //������
		br = bi = 0;
		for (j1 = 0; j1<I; j1++) {      //�������
			br = br + fr[i1*L + j1] * xr[j1] - fi[i1*L + j1] * xi[j1];
			bi = bi + fi[i1*L + j1] * xr[j1] + fr[i1*L + j1] * xi[j1];
		}
		//    DTmat<<br<<"+i"<<bi<<"    "<<fr[i1*L+I]<<"+i"<<fi[i1*L+I]<<"\n";
	}
	//    cout<<"\n";    getch();
	//DTmat<<"�������� � ����������� �������\n";
	for (i1 = 0; i1<I; i1++) {        //������
		br = bi = 0;
		for (j1 = 0; j1<I; j1++) {      //�������
			br = br + ftr[i1*L + j1] * xr[j1] - fti[i1*L + j1] * xi[j1];
			bi = bi + fti[i1*L + j1] * xr[j1] + ftr[i1*L + j1] * xi[j1];
		}
		//    DTmat<<br<<"+i"<<bi<<"    "<<ftr[i1*L+I]<<"+i"<<fti[i1*L+I]<<"\n";
	}
}  //�o��� GAUSSccc
   //*************************************************************************

void GAUSSc3(int I, double ar[], double ai[], double xr[], double xi[]) {
	//GAUSSc3: o���������� �o���� �o�������o� ������� �������� ���������
	//���������: GAUSSc(I,Ar,Ai,Xr,Xi)
	//     const  I          - �������o��� �������
	//     double �r[I*I+I]  -����������� ������� ������������o� �����
	//                        �o���������o�
	//     double �i[I*I+I]  -����������� ������� ����o� ����� �o���������o�
	//     double �r[I]      -����o� ������������o� ����� �o���� �������
	//     double Xi[I]      -����o� ����o� ����� �o���� �������
	//cout<<"I="<<I<<"\n";                     getch();
	int i, j, m, p, q, i1, j1, L;
	double br, bi, cr, ci, dr, di, e;
	L = I + 1;
	for (m = 1; m<I; m++) {
		for (i = m; i<I; i++) {
			p = (m - 1)*I + m + m - 2;
			cr = ar[p];
			ci = ai[p];
			q = i*L + m - 1;
			dr = ar[q];
			di = ai[q];
			e = 1 / (dr*dr + di*di);
			br = (cr*dr + ci*di)*e;
			bi = (ci*dr - cr*di)*e;
			// b=a[(m-1)*I+m-2+m]/a[i*L+m-1];
			for (j = m - 1; j<L; j++) {
				p = i*L + j;
				cr = ar[p];
				ci = ai[p];
				q = (m - 1)*L + j;
				ar[p] = cr*br - ci*bi - ar[q];
				ai[p] = ci*br + cr*bi - ai[q];
				//a[i*L+j]=a[i*L+j]*b-a[(m-1)*L+j];
				// for(i1=0;i1<I;i1++)                                       //���o�
				//  {for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1]<<"  ";cout<<"\n";} //��o�����o����
				//  cout<<"\n";getch();                                       //������
			}
		}
	}  // cout<<"������� �����o����� �������\n"; getch();
	for (j = I - 1; j>0; j--) {
		for (i = j - 1;i >= 0;i--) {
			p = j*L + I;
			cr = ar[p];
			ci = ai[p];
			q = j*L + j;
			dr = ar[q];
			di = ai[q];
			e = 1 / (dr*dr + di*di);
			br = (cr*dr + ci*di)*e;
			bi = (ci*dr - cr*di)*e;
			//b=a[j*L+I]/a[j*L+j];
			p = i*L + I;
			cr = ar[p];
			ci = ai[p];
			q = i*L + j;
			dr = ar[q];
			di = ai[q];
			ar[p] = ar[p] - br*dr + di*bi;
			ai[p] = ai[p] - bi*dr - br*di;
			// a[i*L+I]=a[i*L+I]-b*a[i*L+j];
			//    cout<<"i="<<i<<"  j="<<j<<"\n";    for(i1=0;i1<I;i1++)     //���o�
			//  { for(j1=0;j1<=I;j1++)cout<<a[i1*(I+1)+j1]<<"  ";cout<<"\n";}  //��������o����
			//  cout<<"\n";getch();                                       //������
		}
	}
	//cout<<"�o���:";
	for (i = 0; i<I; i++) {
		p = i*L + I;
		cr = ar[p];
		ci = ai[p];
		q = i*L + i;
		dr = ar[q];
		di = ai[q];
		e = 1 / (dr*dr + di*di);
		xr[i] = (cr*dr + ci*di)*e;
		xi[i] = (ci*dr - cr*di)*e;
		//  cout<< xr[i]<<"  "<<xi[i]<<"\n"; getch();
	}
	//    cout<<"\n";    getch();
}  //�o��� GAUSSc3


   //*************************************************************************
void TREGAUSSC(int const I, long double ar[], long double ai[]) {
	//TREGAUSSC: ���������� ����������� ������� � �o���������� �������������� � ������������ ����
	//���������: TREGAUSSC(I,Ar,Ai)
	//     const  I          - �������o��� �������
	//long double �r[I*I+I]  -����������� ������� ������������o� �����
	//                        �o���������o�
	//long double �i[I*I+I]  -����������� ������� ����o� ����� �o���������o�
	//cout<<"I="<<I<<"\n";                     getch();
	int i, j, m, p, q, i1, j1, L;
	long double br, bi, cr, ci, dr, di, e;
	L = I + 1;
	for (m = 1; m<I; m++) {
		for (i = m; i<I; i++) {
			p = (m - 1)*I + m + m - 2;
			cr = ar[p];
			ci = ai[p];
			q = i*L + m - 1;
			dr = ar[q];
			di = ai[q];
			e = 1 / (dr*dr + di*di);
			br = (cr*dr + ci*di)*e;
			bi = (ci*dr - cr*di)*e;
			// b=a[(m-1)*I+m-2+m]/a[i*L+m-1];
			for (j = m - 1;j<L;j++) {
				p = i*L + j;
				cr = ar[p];
				ci = ai[p];
				q = (m - 1)*L + j;
				ar[p] = cr*br - ci*bi - ar[q];
				ai[p] = ci*br + cr*bi - ai[q];
				//a[i*L+j]=a[i*L+j]*b-a[(m-1)*L+j];
				// for(i1=0;i1<I;i1++)                                       //���o�
				//  {for(j1=0;j1<=I;j1++)cout<<ar[i1*L+j1]<<"  ";cout<<"\n";} //��o�����o����
				//  cout<<"\n";getch();                                       //������
			}
		}
	}  // cout<<"������� �����o����� �������\n"; getch();
} // ����� TREGAUSSC



  //************************************************************************
void GAUSSD(int const I, int const M, float x[], double a[]) {
	//GAUSS: ����������� ������ ������� �������� ��������� � ��������������� ��������������
	//���������: GAUSS(I,A,X)
	//     const  I     - ����������� �������
	//     const  M
	//     double A[]   -����������� ������� �������������
	//            �-      ����� �������

	cout << "I=" << I << "  M=" << M << "\n";
	getch();

	int i, j, m, i1, j1, L;
	double b;
	L = I + 1;
	for (m = 1; m<I; m++) {
		for (i = m; i<I; i++) {
			b = a[(m - 1)*I + m - 2 + m] / a[i*L + m - 1];
			for (j = m - 1; j<L; j++) {
				a[i*L + j] = a[i*L + j] * b - a[(m - 1)*L + j];
				// for(i1=0;i1<I;i1++)                                       //�����
				//  {for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1]<<"  ";cout<<"\n";} //�������������
				//  cout<<"\n";getch();                                       //������
			}
		}
	}

	cout << "������� ����������� �������\n";
	getch();

	for (j = I - 1;j>0;j--) {
		for (i = j - 1;i >= 0;i--) {
			b = a[j*L + I] / a[j*L + j];
			a[i*L + I] = a[i*L + I] - b*a[i*L + j];
			//    cout<<"i="<<i<<"  j="<<j<<"\n";    for(i1=0;i1<I;i1++)     //�����
			//  { for(j1=0;j1<=I;j1++)cout<<a[i1*(I+1)+j1]<<"  ";cout<<"\n";}  //�������������
			//  cout<<"\n";getch();                                       //������
		}
	}

	cout << "�����:";

	for (i = 0; i<I; i++) {
		x[i] = a[i*L + I] / a[i*L + i];
		cout << x[i] << "  ";
		getch();
	}

	cout << "\n";
	getch();
}  //����� GAUSSD


   //*************************************************************************
void TREGAUSSD(int const I, int const M, double a[])
{//TREGAUSSD: -���������� ����������� ������� � ��������������� �������������� � ������������ ����
 //���������: TREGAUSSD(I,A,X)
 //     const  I     - ����������� �������
 //     const  M
 //     double A[]   -����������� ������� �������������
	cout << "I=" << I << "  M=" << M << "\n";                     getch();
	int i, j, m, i1, j1, L;
	double b;
	L = I + 1;
	for (m = 1; m<I; m++) {
		for (i = m; i<I; i++) {
			b = a[(m - 1)*I + m - 2 + m] / a[i*L + m - 1];
			for (j = m - 1; j<L; j++) {
				a[i*L + j] = a[i*L + j] * b - a[(m - 1)*L + j];
				// for(i1=0;i1<I;i1++)                                       //�����
				//  {for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1]<<"  ";cout<<"\n";} //�������������
				//  cout<<"\n";getch();                                       //������
			}
		}
	}

	cout << "������� ����������� �������\n";
	getch();

}  //����� TREGAUSSD
#endif __GAUSS_LIB

