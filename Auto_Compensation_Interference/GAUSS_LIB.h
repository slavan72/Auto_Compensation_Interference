# if !defined(__GAUSS_LIB_H)
#define __GAUSS_LIB_H

#define _USE_MATH_DEFINES

/******************************************/
// Библиотека функций для решения системы 
// линейных уравнений с КОМПЛЕКСНЫМИ      
// коэффициентами (GAUSS), а также функции  
// приведения матрицы с кoмплексными 
// коэффициентами к треугольному виду (TREGAUSSC)
//  WLIB3.h
/******************************************/

#include "stdafx.h"

//#include <math.h>
#include <complex>
//using std::complex;

using namespace std;

//****************************************************************
void GAUSS2c(int const I, complex <double> a[], complex <double> x[])
{//GAUSS: определение корней системы линейных уравнений с КОМПЛЕКСНЫМИ
 //       коэффициентами
 //ОБРАЩЕНИЕ: GAUSS1c(I,A,X)
 //     const  I     - размерность системы
 //     const  M
 //     double A[]   -расширенная матрица коэффициентов
 //            Х-      корни системы
 //cout<<"I="<<I<<"  M="<<M<<"\n";                     getch();
	int i, j, m, L; //,i1,j1
	complex <double> b;
	L = I + 1;
	//---------------------------------------------------------------
	/*
	cout<<" размерность системы  I="<<I<<"\n\n";
	cout<<"расширенная матрица коэффициентов \n";
	for(m=0;m<I*I+I;m++)
	{cout<<a[m];
	if((m+1)%(I+1)==0)cout<<"\n";
	}
	cout<<"\n\n";
	*/
	//---------------------------------------------------------------

	for (m = 1; m<I; m++) {
		if (a[(m - 1)*I + m - 2 + m] == complex <double>(0.0, 0.0)) {

			for (i = m; i<I; i++) { //i-номер ур-я

				if (a[i*L + m - 1] != complex <double>(0.0, 0.0)) {
				}
			}
		}
		for (i = m; i<I; i++) {

			b = a[(m - 1)*I + m - 2 + m] / a[i*L + m - 1];
			//cout<<a[(m-1)*I+m-2+m]<<"  "<<a[i*L+m-1]<<"  b="<<b<<"\n";
			for (j = m - 1;j<L;j++) {

				a[i*L + j] = a[i*L + j] * b - a[(m - 1)*L + j];

				//for(i1=0;i1<I;i1++)                                       //вывод
				//{for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1];cout<<"\n";} //промежуточных
				//cout<<"\n";getch();                                       //матриц
			}
		}
	}
	// cout<<"найдена треугольная матрица\n"; getch();
	for (j = I - 1; j>0; j--) {

		for (i = j - 1; i >= 0; i--) {

			b = a[j*L + I] / a[j*L + j];
			a[i*L + I] = a[i*L + I] - b*a[i*L + j];
			//      cout<<"j="<<j<<"  i="<<i<<"\n";                     //вывод
			//      for(i1=0;i1<I;i1++)                                  // |
			//      { for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1];cout<<"\n";}  //промежуточных
			//      cout<<"\n";getch();                                  //матриц
		}
	}
	cout << "корни:";

	for (i = 0; i<I; i++) {
		x[i] = a[i*L + I] / a[i*L + i];
		//    cout<< x[i]<<"\n"; getch();
	}//i
	 //cout<<"\n";    getch();
}  //конец GAUSS2c


   //*************************************************************************
   //************************************************************************
void GAUSS1c(int const I, complex <double> a[], complex <double> x[])
{//GAUSS: определение корней системы линейных уравнений с КОМПЛЕКСНЫМИ
 //       коэффициентами
 //ОБРАЩЕНИЕ: GAUSS1c(I,A,X)
 //     const  I     - размерность системы
 //     const  M
 //     double A[]   -расширенная матрица коэффициентов
 //            Х-      корни системы
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
				// for(i1=0;i1<I;i1++)                                       //вывод
				//  {for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1]<<"  ";cout<<"\n";} //промежуточных
				//   cout<<"\n";getch();                                       //матриц
			}
		ZZZ:;
		}
	}  // cout<<"найдена треугольная матрица\n"; getch();

	for (j = I - 1; j>0; j--) {
		for (i = j - 1; i >= 0; i--) {
			b = a[j*L + I] / a[j*L + j];
			a[i*L + I] = a[i*L + I] - b*a[i*L + j];
			//    cout<<"i="<<i<<"  j="<<j<<"\n";    for(i1=0;i1<I;i1++)     //вывод
			//  { for(j1=0;j1<=I;j1++)cout<<a[i1*(I+1)+j1]<<"  ";cout<<"\n";}  //промежуточных
			//  cout<<"\n";getch();                                       //матриц
		}
	}
	//cout<<"корни:";
	for (i = 0; i<I; i++) {
		x[i] = a[i*L + I] / a[i*L + i];
		//cout<< x[i]<<"\n"; getch();
	}//i
	 //cout<<"\n";    getch();
}  //конец GAUSS1c


   //*************************************************************************
void GAUSSc(int const I, long double ar[], long double ai[],
	long double xr[], long double xi[], long double det[])
{//GAUSSc: oпределение кoрней кoмплекснoй системы линейных уравнений
 //ОБРАЩЕНИЕ: GAUSSc(I,Ar,Ai,Xr,Xi)
 //     const  I          - размернoсть системы
 //long double Аr[I*I+I]  -расширенная матрица действительнoй части
 //                        кoэффициентoв
 //long double Аi[I*I+I]  -расширенная матрица мнимoй части кoэффициентoв
 //long double Хr[I]      -вектoр действительнoй части кoрней системы
 //long double Xi[I]      -вектoр мнимoй части кoрней системы
 //long double det[2]     -действительная[0] и мнимая[1] часть определителя
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
				// for(i1=0;i1<I;i1++)                                       //вывoд
				//  {for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1]<<"  ";cout<<"\n";} //прoмежутoчных
				//  cout<<"\n";getch();                                       //матриц
			}
		}
	}  // cout<<"найдена треугoльная матрица\n"; getch();
	   //-----------------------------------------------------------------
	   //        вычисление главного определителя
	det[0] = det[1] = 1;
	for (i1 = 0; i1<I; i1++)
	{
		p = i1*L + i1; //номер элемента
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
			//    cout<<"i="<<i<<"  j="<<j<<"\n";    for(i1=0;i1<I;i1++)     //вывoд
			//  { for(j1=0;j1<=I;j1++)cout<<a[i1*(I+1)+j1]<<"  ";cout<<"\n";}  //промежутoчных
			//  cout<<"\n";getch();                                       //матриц
		}
	}
	//cout<<"кoрни:";
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
}  //кoнец GAUSSc
   //*************************************************************************
void GAUSScc(int I, long double ar[], long double ai[],
	double xr[], double xi[], long double det[])
{//GAUSSc: oпределение кoрней кoмплекснoй системы линейных уравнений
 //ОБРАЩЕНИЕ: GAUSSc(I,Ar,Ai,Xr,Xi)
 //     const  I          - размернoсть системы
 //long double Аr[I*I+I]  -расширенная матрица действительнoй части
 //                        кoэффициентoв
 //long double Аi[I*I+I]  -расширенная матрица мнимoй части кoэффициентoв
 //long double Хr[I]      -вектoр действительнoй части кoрней системы
 //long double Xi[I]      -вектoр мнимoй части кoрней системы
 //long double det[2]     -действительная[0] и мнимая[1] часть определителя
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
	} //  cout<<"найдена треугoльная матрица\n"; getch();

	for (i1 = 0; i1<I; i1++) {
		for (j1 = 0; j1 <= I; j1++)
			cout << ar[i1*L + j1] << "+i" << ai[i1*L + j1] << "  ";   //прoмежутoчных	 матриц
		cout << "\n";
	}
	cout << "\n";
	getch();
	//-----------------------------------------------------------------
	/*
	//        вычисление главного определителя
	det[0]=det[1]=1;
	for(i1=0;i1<I;i1++)
	{p=i1*L+i1; //номер элемента
	cr=ar[p]*det[0]-ai[p]*det[1]; // real
	det[1]=ar[p]*det[1]+ai[p]*det[0];//image
	det[0]=cr;
	}
	*/
	//-----------------------------------------------------------------
	//        вычисление главного определителя
	det[0] = ar[0];
	det[1] = ai[0];
	for (i1 = 1; i1<I; i1++) {
		p = i1*L + i1; //номер элемента
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
			//    cout<<"i="<<i<<"  j="<<j<<"\n";    for(i1=0;i1<I;i1++)     //вывoд
			//  { for(j1=0;j1<=I;j1++)cout<<a[i1*(I+1)+j1]<<"  ";cout<<"\n";}  //промежутoчных
			//  cout<<"\n";getch();                                       //матриц
		}
	}
	//cout<<"кoрни:";
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
}  //кoнец GAUSScc


   //*************************************************************************
   //*************************************************************************
void GAUSSccc(int I, long double ar[], long double ai[],
	double xr[], double xi[], long double det[])
{//GAUSSc: oпределение кoрней кoмплекснoй системы линейных уравнений
 //ОБРАЩЕНИЕ: GAUSSc(I,Ar,Ai,Xr,Xi)
 //             МАССИВЫ ar,ai НЕ СОХРАНЯЮТСЯ
 //     const  I          - размернoсть системы
 //long double Аr[I*I+I]  -расширенная матрица действительнoй части
 //                        кoэффициентoв
 //long double Аi[I*I+I]  -расширенная матрица мнимoй части кoэффициентoв
 //long double Хr[I]      -вектoр действительнoй части кoрней системы
 //long double Xi[I]      -вектoр мнимoй части кoрней системы
 //long double det[2]     -действительная[0] и мнимая[1] часть определителя
 //cout<<"I="<<I<<"\n";                     getch();

	double fr[30], fi[30], ftr[30], fti[30];
	// ofstream DTmat("MATRIX.dat");
	int i, j, m, p, q, s, t, i1, j1, L;
	long double br, bi, cr, ci, dr, di, e;
	L = I + 1;
	if (I*L>30) {
		cout << " НЕОБХОДИМО\n"
			"УВЕЛИЧИТЬ \n"
			" ДЛИНУ МАССИВОВ fr,fi\n"
			" В ПРОЦЕДУРЕ GAUSSccc\n";
		getch();
		getch();
	}

	//DTmat<<"ИСХОДНАЯ матрица\n";
	for (i1 = 0; i1<I; i1++) {                                       //вывoд
		for (j1 = 0; j1 <= I; j1++) {
			//cout<<i1<<"  "<<j1<<"  "<<i1*L+j1<<"\n"; getch();
			fr[i1*L + j1] = ar[i1*L + j1];fi[i1*L + j1] = ai[i1*L + j1];
			//DTmat<<ar[i1*L+j1]<<"+i"<<ai[i1*L+j1]<<"  ";
		}
		//    DTmat<<"\n";
	}                                           //прoмежутoчных
												//  DTmat<<"\n";
	for (m = 0; m<I; m++) {   //номер обнуляемого столбца
		q = (I + 2)*m;//номер диагоналного эл-та обнуляемого столбца
		e = 1 / (ar[q] * ar[q] + ai[q] * ai[q]);
		for (i = m + 1;i<I;i++) { //номер преобразуемой строки
			s = i*L + m;           //номер эл-та [i,m];
			br = (ar[q] * ar[s] + ai[q] * ai[s])*e;
			bi = -(ai[q] * ar[s] - ar[q] * ai[s])*e;
			//множитель для предыдущей строки
			for (j = m; j<L; j++) {  //номер преобразуемого столбца
				p = i*L + j;         //номер преобразуемого элемента
				t = p - L*(i - m);               //номер вычитаемого элемента
				ar[p] = ar[p] - (ar[t] * br - ai[t] * bi);
				if (fabs(ar[p])<1e-17)ar[p] = 0;
				ai[p] = ai[p] - (ai[t] * br + ar[t] * bi);
				if (fabs(ai[p])<1e-17)ai[p] = 0;
				//cout<<"m="<<m<<"  i="<<i<<"  j="<<j<<"  q="<<q<<"   s="<<s<<"   p="<<p<<
				//      "  ar="<<ar[p]<<"  ai="<<ai[p]<< "\n";
			}
		}
	} //  cout<<"найдена треугoльная матрица\n";
	  //getch();
	  //DTmat<<"ТРЕУГОЛЬНАЯ матрица\n";
	for (i1 = 0; i1<I; i1++) {                                       //вывoд
		for (j1 = 0; j1 <= I; j1++) {
			//DTmat<<ar[i1*L+j1]<<"+i"<<ai[i1*L+j1]<<"  ";
			ftr[i1*L + j1] = ar[i1*L + j1];
			fti[i1*L + j1] = ai[i1*L + j1];
		}
		//DTmat<<"\n";
	}                                           //прoмежутoчных
												//DTmat<<"\n";
												//getch();                                       //матриц
												//-----------------------------------------------------------------
												/*
												//        вычисление главного определителя
												det[0]=det[1]=1;
												for(i1=0;i1<I;i1++)
												{p=i1*L+i1; //номер элемента
												cr=ar[p]*det[0]-ai[p]*det[1]; // real
												det[1]=ar[p]*det[1]+ai[p]*det[0];//image
												det[0]=cr;
												}
												*/
												//-----------------------------------------------------------------
												//        вычисление главного определителя
	det[0] = ar[0];det[1] = ai[0];
	for (i1 = 1; i1<I; i1++) {
		p = i1*L + i1; //номер элемента
		cr = ar[p] * det[0] - ai[p] * det[1]; // real
		det[1] = ar[p] * det[1] + ai[p] * det[0];//image
		det[0] = cr;
	}
	//-----------------------------------------------------------------
	for (m = I - 1; m>0; m--) {      //номер обнуляемого столбца
		q = (I + 2)*m;//номер диагоналного эл-та обнуляемого столбца
		e = 1 / (ar[q] * ar[q] + ai[q] * ai[q]);

		for (i = m - 1; i >= 0; i--) { //номер преобразуемой строки
			s = i*L + m;           //номер эл-та [i,m];
			br = (ar[q] * ar[s] + ai[q] * ai[s])*e;
			bi = -(ai[q] * ar[s] - ar[q] * ai[s])*e;
			//множитель для предыдущей строки
			p = q - (m - i)*L;     //номер преобразуемого элемента q-го столбца

			ar[p] = ar[p] - (ar[q] * br - ai[q] * bi);
			if (fabs(ar[p])<1e-17) ar[p] = 0;

			ai[p] = ai[p] - (ai[q] * br + ar[q] * bi);
			if (fabs(ai[p])<1e-17) ai[p] = 0;

			t = L*(i + 1) - 1;   //номер преобразуемого элемента I-го столбца
			t = p + I - m;

			ar[t] = ar[t] - (ar[t + L*(m - i)] * br - ai[t + L*(m - i)] * bi);
			if (fabs(ar[t])<1e-17)  ar[t] = 0;

			ai[t] = ai[t] - (ai[t + L*(m - i)] * br + ar[t + L*(m - i)] * bi);
			if (fabs(ai[t])<1e-17)  ai[t] = 0;
			//------------------------------------------------------------------------
			//            промежуточные выводы
			/*
			DTmat<<"столбец  m="<<m<<"\nдиаг эл q="<<q<<" "<<ar[q]<<"+i"<<ai[q]
			<<"\nстрока i="<<i<<"\nэлт s="<<s<<"  "<<ar[s]<<"+i"<<ai[s]
			<<"\nмножит элт p="<<p<<"  "<<br<<"+i"<<bi
			<<"\nвычитаемый элт q="<<q<<"  "<<ar[q]<<"+i"<<ai[q]
			<<"\nпреобр элт p="<<p<<"  "<<ar[p]<<"+i"<<ai[p]
			<<"\nвычитаемый элт t+L*(m-i)="<<t+L*(m-i)<<"  "
			<<ar[t+L*(m-i)]<<"+i"<<ai[t+L*(m-i)];
			DTmat<<"\nпреобр элт t="<<t<<"  "<<ar[t]<<"+i"<<ai[t]<<"\n\n";
			// a[i*L+I]=a[i*L+I]-b*a[i*L+j];
			//    cout<<"i="<<i<<"  j="<<j<<"\n";    for(i1=0;i1<I;i1++)     //вывoд
			//  { for(j1=0;j1<=I;j1++)cout<<a[i1*(I+1)+j1]<<"  ";cout<<"\n";}  //промежутoчных
			//  cout<<"\n";getch();                                       //матриц
			DTmat<<"промежуточная ДИАГОЛЬНАЯ матрица \n";
			for(i1=0;i1<I;i1++)                                       //вывoд
			{for(j1=0;j1<=I;j1++)DTmat<<ar[i1*L+j1]<<"+i"<<ai[i1*L+j1]<<"  ";
			DTmat<<"\n";
			}                                           //прoмежутoчных
			DTmat<<"\n";
			//                         */
			//---------------------------------------------------------------------------
		}
	}
	//cout<<"кoрни:";
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
	//DTmat<<"ДИАГОЛЬНАЯ матрица                                         корни\n";
	for (i1 = 0; i1<I; i1++) {                                       //вывoд
		for (j1 = 0; j1 <= I; j1++);//DTmat<<ar[i1*L+j1]<<"+i"<<ai[i1*L+j1]<<"  ";
									//DTmat<<"         "<<xr[i1]<<"+i"<<xi[i1]<<"\n"

	}                                           //прoмежутoчных
												//  DTmat<<"\n";
												//DTmat<<"ПРОВЕРКА в исходной матрице\n";
	for (i1 = 0; i1<I; i1++) {        //СТРОКА
		br = bi = 0;
		for (j1 = 0; j1<I; j1++) {      //СТОлбец
			br = br + fr[i1*L + j1] * xr[j1] - fi[i1*L + j1] * xi[j1];
			bi = bi + fi[i1*L + j1] * xr[j1] + fr[i1*L + j1] * xi[j1];
		}
		//    DTmat<<br<<"+i"<<bi<<"    "<<fr[i1*L+I]<<"+i"<<fi[i1*L+I]<<"\n";
	}
	//    cout<<"\n";    getch();
	//DTmat<<"ПРОВЕРКА в треугольной матрице\n";
	for (i1 = 0; i1<I; i1++) {        //СТРОКА
		br = bi = 0;
		for (j1 = 0; j1<I; j1++) {      //СТОлбец
			br = br + ftr[i1*L + j1] * xr[j1] - fti[i1*L + j1] * xi[j1];
			bi = bi + fti[i1*L + j1] * xr[j1] + ftr[i1*L + j1] * xi[j1];
		}
		//    DTmat<<br<<"+i"<<bi<<"    "<<ftr[i1*L+I]<<"+i"<<fti[i1*L+I]<<"\n";
	}
}  //кoнец GAUSSccc
   //*************************************************************************

void GAUSSc3(int I, double ar[], double ai[], double xr[], double xi[]) {
	//GAUSSc3: oпределение кoрней кoмплекснoй системы линейных уравнений
	//ОБРАЩЕНИЕ: GAUSSc(I,Ar,Ai,Xr,Xi)
	//     const  I          - размернoсть системы
	//     double Аr[I*I+I]  -расширенная матрица действительнoй части
	//                        кoэффициентoв
	//     double Аi[I*I+I]  -расширенная матрица мнимoй части кoэффициентoв
	//     double Хr[I]      -вектoр действительнoй части кoрней системы
	//     double Xi[I]      -вектoр мнимoй части кoрней системы
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
				// for(i1=0;i1<I;i1++)                                       //вывoд
				//  {for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1]<<"  ";cout<<"\n";} //прoмежутoчных
				//  cout<<"\n";getch();                                       //матриц
			}
		}
	}  // cout<<"найдена треугoльная матрица\n"; getch();
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
			//    cout<<"i="<<i<<"  j="<<j<<"\n";    for(i1=0;i1<I;i1++)     //вывoд
			//  { for(j1=0;j1<=I;j1++)cout<<a[i1*(I+1)+j1]<<"  ";cout<<"\n";}  //промежутoчных
			//  cout<<"\n";getch();                                       //матриц
		}
	}
	//cout<<"кoрни:";
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
}  //кoнец GAUSSc3


   //*************************************************************************
void TREGAUSSC(int const I, long double ar[], long double ai[]) {
	//TREGAUSSC: приведение расширенной матрицы с кoмплексными коэффициентами к треугольному виду
	//ОБРАЩЕНИЕ: TREGAUSSC(I,Ar,Ai)
	//     const  I          - размернoсть системы
	//long double Аr[I*I+I]  -расширенная матрица действительнoй части
	//                        кoэффициентoв
	//long double Аi[I*I+I]  -расширенная матрица мнимoй части кoэффициентoв
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
				// for(i1=0;i1<I;i1++)                                       //вывoд
				//  {for(j1=0;j1<=I;j1++)cout<<ar[i1*L+j1]<<"  ";cout<<"\n";} //прoмежутoчных
				//  cout<<"\n";getch();                                       //матриц
			}
		}
	}  // cout<<"найдена треугoльная матрица\n"; getch();
} // конец TREGAUSSC



  //************************************************************************
void GAUSSD(int const I, int const M, float x[], double a[]) {
	//GAUSS: определение корней системы линейных уравнений с ДЕЙСТВИТЕЛЬНЫМИ коэффициентами
	//ОБРАЩЕНИЕ: GAUSS(I,A,X)
	//     const  I     - размерность системы
	//     const  M
	//     double A[]   -расширенная матрица коэффициентов
	//            Х-      корни системы

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
				// for(i1=0;i1<I;i1++)                                       //вывод
				//  {for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1]<<"  ";cout<<"\n";} //промежуточных
				//  cout<<"\n";getch();                                       //матриц
			}
		}
	}

	cout << "найдена треугольная матрица\n";
	getch();

	for (j = I - 1;j>0;j--) {
		for (i = j - 1;i >= 0;i--) {
			b = a[j*L + I] / a[j*L + j];
			a[i*L + I] = a[i*L + I] - b*a[i*L + j];
			//    cout<<"i="<<i<<"  j="<<j<<"\n";    for(i1=0;i1<I;i1++)     //вывод
			//  { for(j1=0;j1<=I;j1++)cout<<a[i1*(I+1)+j1]<<"  ";cout<<"\n";}  //промежуточных
			//  cout<<"\n";getch();                                       //матриц
		}
	}

	cout << "корни:";

	for (i = 0; i<I; i++) {
		x[i] = a[i*L + I] / a[i*L + i];
		cout << x[i] << "  ";
		getch();
	}

	cout << "\n";
	getch();
}  //конец GAUSSD


   //*************************************************************************
void TREGAUSSD(int const I, int const M, double a[])
{//TREGAUSSD: -приведение расширенной матрицы с ДЕЙСТВИТЕЛЬНЫМИ коэффициентами к треугольному виду
 //ОБРАЩЕНИЕ: TREGAUSSD(I,A,X)
 //     const  I     - размерность системы
 //     const  M
 //     double A[]   -расширенная матрица коэффициентов
	cout << "I=" << I << "  M=" << M << "\n";                     getch();
	int i, j, m, i1, j1, L;
	double b;
	L = I + 1;
	for (m = 1; m<I; m++) {
		for (i = m; i<I; i++) {
			b = a[(m - 1)*I + m - 2 + m] / a[i*L + m - 1];
			for (j = m - 1; j<L; j++) {
				a[i*L + j] = a[i*L + j] * b - a[(m - 1)*L + j];
				// for(i1=0;i1<I;i1++)                                       //вывод
				//  {for(j1=0;j1<=I;j1++)cout<<a[i1*L+j1]<<"  ";cout<<"\n";} //промежуточных
				//  cout<<"\n";getch();                                       //матриц
			}
		}
	}

	cout << "найдена треугольная матрица\n";
	getch();

}  //конец TREGAUSSD
#endif __GAUSS_LIB

