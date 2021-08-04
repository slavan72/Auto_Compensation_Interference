// Auto_Compensation_Interference.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"
//#include "Console_String_LIB.h"
#include "Simple_math_LIB.h"
//#include "GAUSS_LIB.h"
//#include "WLIB4.h"
#include "F4D_L.H"
#include "MyGraphics.h"

#include <iostream>
#include <fstream>

int const Nr = 3;//601;//81;//количество реализаций начиная с первой
				 //int const Nc = 3;//количество каналов  (включая основной)
int const NC = 4;//максимальное кол-во каналов

int tstep = 125;//208;  //длительность шага  нс
int Np = 2;//количество помех
int Nc = 3;//количество каналов  (включая основной)

double tzam = 1. / (fmax0 - fmin0); //период замеров для colibrovca
double zader[5][40];//, divider;
float lam = 300 / 9.5;//длина волны[mм]
float ashum = 50e-6; //напряжение собст шума(вольт)

double af_1, uv_f, u_f;  //int n_f;//,m_f; //используется внутри функций

int const Ntipfv = 2;
int const Nfv = 15;
float Q;   //


float pomeha[4][3] =//номер помехи,коорд-ты в АЗУМ и мощность помехи
{
	1.0, -1.33, 23.5,//48.840,//113.94,//54.463,66.716, //ум,аз,P  1-oй помехи
	1.0, 1.33, 55.5,//48.840,//43.786,
	1.0, 31.9, 54.7,//43.786,
	1.0, -14., 42.716
};

float azumant[5][3] = //номер антенны,коорд луча в АЗУМ и тип решет
{ 1., .0, 0, //UM,AZ,тип решет  канала[0](основного)
1., -1.33, 2,//-17.5,            //канал[1]
1., 1.33, 1,//-17.5,            //канал[1]
1.,31.9,1,//-17.5,            //канал[1]
1.,31.9,1,//-17.5,            //канал[3]
		  //  1.,6.3,1,//-17.5,            //канал[2]
		  //  1.,2.,1//-17.5            //канал[4]
};

double dt = 1e-9;//*125./126.;
double uints[];
double uintc[];

struct FVTIP {
	double ga[Ntipfv];
	double gb[Ntipfv];
	double a[Ntipfv];
	double b[Ntipfv];
	float k[Ntipfv];
}FVST, FVCT;

struct FVNOM {
	double u0[Nfv];
	double ua[Nfv];
	double ub[Nfv];
	double du[Nfv];
	double un[Nfv];
}FVSN, FVCN;

struct FVp {
	float f0[Ntipfv];
	float as[Ntipfv];
	float bs[Ntipfv];
	float ac[Ntipfv];
	float bc[Ntipfv];
	float h[Ntipfv];
}FVpsi;

//****************************************************
//                   КОНТУР
void QFO(float f1, float f2, int n, float& fo, float& Q);

void FVHT(float dt, float f0, int nomfvtip, FVTIP& FVST, FVTIP& FVCT, FVp& FVpsi);//инициализируем структуры FVST, FVCT, FVpsi

void UFVn(FVNOM& FVSN, FVNOM& FVCN); // обнуляем структуры FVSN, FVCN

void KHT(float dt, float f0, float Q, float W, float Ku, int nomkontip, KONTIP& KT);//инициализируем структуру KT

void DIFCHT(float dt, float tau, int tip, DIFCTIP& DIFCT);	//инициализируем структуру DIFCT

void UFVn(FVNOM& FVSN, FVNOM& FVCN);// обнуляем структуры FVSN, FVCN

void UFV(int nom, int tip, FVTIP& FVST, FVTIP&FVCT, FVNOM& FVSN, FVNOM& FVCN, double uv, double& us, double& uc);//инициализируем структуры

void UDIFCn(DIFCNOM& DIFCN);

void ANTCAN(float naklon[2], float geom[5][2], float azumant[5][3], float pomeha[4][3], int Np, int Nc, double zad[4][5], double pole[4][5][2]);

void SIGNAL(double dt, double fo, long int n, int Nsl, double kalibr, int Nc, int Np, double v[4][5][2], double af[4][5][2],
	double tns[5], double vs[5][2], double afs[5][2],
	double zad[4][5], double pole[4][5][2],
	double usignal[5]);

void UPOMCAN(double dt, double fo, long int n, int Nsl, double zad, double v[2], double af[2], double pole[2], double&u);

void SDVIG(int Nc, double usignal[], double usignals[], double usignalc[]);

void ZADERG(int nom, int nzad, long int n, double u1, double& u2);

void SHUM(double dt, double fo, long int n, int Nsl, double tns, double vs[2], double afs[2], double& u);

void REGUS(int L, int Nc, double usignals[], double usignalc[], double uints[], double uintc[], double uws[], double uwc[], double Wc[], double Ws[]);

void PROCESSORdig1(int Nc, int Np, int nstep, double udigs[], double udigc[], double ws[], double wc[], double Wds[], double Wdc[], double Pnorm);

void ANALCAN(int Nc, double usignal[], double udig[]);

void ANALCANos(double usigm, double& usigm1);


int main()
{
	setlocale(LC_ALL, "Russian");
	/***COMMENT***/
	float X[100], Y[100];//,Z[600];//,Z1[600];//,Z2[600];
	double pc[NC], pcc[NC], pcr[NC], pci[NC], kor[NC];//
	double  ucdigzad[5][25], usdigzad[5][25];
	long int L, q;//
				  //int Nusred = 1;//количество усредняемых замеров
	int n, m, nc, nstep, ntakt;//
	float Ksl[5];
	double ksl;
	float zagr[Nr], zagr1[Nr];
	double t;//, tkon;//tkon-момент окончания адаптации
	int Nsl, Nn, Nk, Ngraf;	//Nust, Nw,
	int Nadap = 0;
	double tzam = 1. / (fmax0 - fmin0); //период замеров для colibrovca
										/***COMMENT***/

										//double foo=6e7; //средняя частота сигнала на ПЧ

										//		       /*
										//             Длительность адаптации   70 mks
										//double tw=70.9e-6;  //+длительность нарастания W


	double tadap = 70e-6; //+длительность адаптации
						  //int Nstep = (int)(tadap*1e9/tstep);//576;
	double tn = 75e-6;//+начало усреднения сигнала
	double tgraf = tn + 10e-6;//130e-6;//длительность отображения на графиках
	double tk = tn + 50e-6;//+окончание усреднения сигнала

						   /***COMMENT***/
	double Wc[NC - 1], Wdc[NC - 1], Ws[NC - 1], Wds[NC - 1], uws[NC - 1], uwc[NC - 1], uints[NC - 1], uintc[NC - 1], usdig[NC], ucdig[NC];//,
	double	uints_zad[NC - 1], uintc_zad[NC - 1]; //	wsigc[NC-1],wsigs[NC-1];

												  //---------------------------------------------------------------------
	float naklon[2] = { 30.,0. };//наклон и поворот полотна град
	float  geom[NC][2] =      //антенна,коорд,
	{
		0.,0.,   //фазовый центр канала[0](основного)
				 //-35.81,-51.01,
				 //	0.,-51.01,
		0.,-51.01,
		-29.4,-51.01,   //канал[1](компенс)  x=-35.81
		31.06,-51.01,
	};
	//******************************************************************
	/***COMMENT***/

	ofstream DTdiag("DIAG.dat", ios::app);
	ofstream DTw_dig("WES_dig.dat");
	ofstream DTz("ZAG.dat");
	ofstream DTze("ZAGend.dat");
	ofstream DTwe("WESend.dat");
	ofstream DTzs("ZAGsre.dat");
	ofstream DTp("PROPIS.dat");
	ofstream DTkk("koefkor.dat");
	ofstream DTzsa("ZAGsred.dar", ios::app);
	//ofstream DTzsa("ZAG_k04.dat",ios::app);
	ofstream DTwsa("WESsred.dar", ios::app);
	ofstream DTmem("memor.dat");
	ofstream DTukor("Ucor.dat", ios::app);
	// ofstream Dpzu("Pzu1.dat");
	ofstream DTw("WES.dat");

	cout << "  I_42 = " << (int)I_42 << "  J_42[0] = " << (int)J_42[0] << "\n";//char I_42=16 - строки; char J_42[16] - массив; F4D_L.h
	_getch();
	/***COMMENT***/
	//int px1, px2, py1, py2;  
	//float x1, x2, y1, y2, dx, dy;
	double pnorm, pn, pk, pk1, pn1, pnorm1, pn2, podav, podav1, podav2, podav21, podav3;
	int zadacp, jjj;
	double a;

	double usignal[5], usignalc[5], usignals[5], us[7];//,uc[7], ufds[5], ufdc[5];
	double usigm, usigm1, usigm2, usigm3;//, usignalzad,  uzad;//,ur[5],ui[5],CF[5],SF[5];
										 /***COMMENT***/

										 //cout<<" dx="<<dx<<"\n"; getch();getch();
	cout << "TIPFAR=" << TIPFAR << "  G=" << G[2];//int const TIPFAR=3; - количество типов реш	 F4D_L.h
	cout << "\n";

	// ifstream READ;
	// READ.open("WES.dat");
	// ofstream TX("WWW6.txt");
	// ofstream DT1("WWW1.dat");
	// ofstream DT2("WWW2.dat",ios::app); //-дозапись в конец файла
	// ofstream DT2("WWW2.dat");      // -cтирание старой информ. и запись новой
	//**********************************************************************
	/*
	for(n=0;n<480;n++)
	{m=(int)(floor(.044*exp(-(n)/500.)*256.+.5));
	cout<<n<<" "<<m<<"\n";    //1 канал  1 помеха
	if(n==0)DTmem<<"0000";
	if(n<10&&n!=0)DTmem<<"000"<<n;
	if(n<100&&n>9)DTmem<<"00"<<n;
	if(n>=100)DTmem<<"0"<<n;
	DTmem<<":";
	if(m==0)DTmem<<"0000";
	if(m<10&&m!=0)DTmem<<"000"<<m;
	if(m==10)DTmem<<"000A";
	if(m>=11)DTmem<<"000B";
	DTmem<<";\n";
	getch();
	}
	*/
	//**********************************************************************
	//unsigned int Nelem_kor[TIPFAR];
	// NELEM_kor(J_kor0,I_kor0,Nelem_kor[0]);amp[0]=(unsigned char*)calloc(Nelem_kor[0],1);
	// NELEM_kor(J_kor1,I_kor1,Nelem_kor[1]);amp[1]=(unsigned char*)calloc(Nelem_kor[1],1);
	// NELEM_kor(J_kor2,I_kor2,Nelem_kor[2]);amp[2]=(unsigned char*)calloc(Nelem_kor[2],1);

	unsigned int Nelem_4[TIPFAR];

	NELEM_F4(I_40, J_40, Nelem_4[0]);	amp[0] = (unsigned char*)calloc(Nelem_4[0], 1);  //cout<<"calloc Nelem_4[0] =" <<amp[0]<<"\n";//отладка
	NELEM_F4(I_41, J_41, Nelem_4[1]);	amp[1] = (unsigned char*)calloc(Nelem_4[1], 1);  //cout<<"calloc Nelem_4[1] =" <<amp[1]<<"\n";//отладка
	NELEM_F4(I_42, J_42, Nelem_4[2]);	amp[2] = (unsigned char*)calloc(Nelem_4[2], 1);  //cout<<"calloc Nelem_4[2] =" <<amp[2] <<"\n";//отладка

	ampn[0] = 0; //int ampn[TIPFAR] -  F4D_L.h

	for (n = 1; n < TIPFAR; n++)	ampn[n] = ampn[n - 1] + ampn[n];// 0+0=0 ???	cout<<"ampn[n] =" <<ampn[n] <<"\n";//отладка
	for (n = 0; n < TIPFAR; n++)	cout << "  Nelem_4[" << n << "] = " << Nelem_4[n] << "\n";// здесь выводятся значения полученные в  NELEM_F4

	_getch();

	for (n = 0; n < TIPFAR; n++) AMP400(n, amp, normir);

	//cout<<DISKR(-52.9,32.,2.)<<"\n";  getch();

	//**************************************************************************
	//double a,b,c;
	/***COMMENT***/
	int nr, j, i, np;//,pp;	
	double kalibr, pnormS, Pnorm;//,u1, u2, u3, uv, uv1;//,us,uc,uv;
								 //double res, re0, re1, ims, im0, im1, ww;     
								 //float diag1, diag2, diag3;
								 //double diag3c, diag3s, psiww, wwc, wws;
								 //L=0;  
								 //long int Lgr; 
								 //int Npoint;
								 /***COMMENT***/

	DTdiag << "1234567\n";
	DTze << Nr << "\n";
	DTzs << Nr << "\n";
	DTzsa << Nr << "\n";
	DTp << Nr << "\n";
	DTwe << Nr << "  " << Nc << "\n";
	DTwsa << Nr << "  " << Nc << "\n";
	DTz << Nr << "\n";
	DTw << Nr << "\n";
	DTw_dig << Nr << "\n";
	DTkk << Nr << "\n";

	//int tipfar[NC]={0,1,2,1,1};
	//1321,1132,1123,1312,   2213,2231, 3312,3321, 1213,1231,
	//2131,2113,  2311,  2312,2321,
	//3112,3121   3221,3212    3312,3321
	//**********************************************************************

	//*********************************************************************
	/*
	float const X[Q]={-23.24,-10.9,11.62,23.24};       //х-кooрдинаты     фазoвых центрoв
	float const Y[Q]={-32.65,-32.65,-32.65,-32.65};    //y-кooрдинаты    кoмпенсаторoв(ламбда)
	*/
	//********************************************************************
	//	   /*

	//int col[NC]={11,12,9,10,13};  int colo;
	//int Iy[TIPFAR]={80,16,16};//kol-vo strok
	//int Jx[TIPFAR]={96,8,8};//kol-vo stolbcov

	/***COMMENT***/
	//int colo;
	//int col[NC]={11,12,9,10};  
	//int Iy[TIPFAR]={112,16,16};//kol-vo strok
	//int Jx[TIPFAR]={96,16,48};//kol-vo stolbcov
	/***COMMENT***/

	//Jx[0]=96;Jx[1]=16;Jx[2]=48;//Jx42;

	//cout<<"Jx = "<<Jx[0]<<"\n";
	//			 /*
	//********************************************************************
	/***COMMENT***/
	double v[4][5][2]; //помеха,канал,компонента
	double af[4][5][2]; //помеха,канал,амплитуда и фаза
	double zad[4][5];   //помеха,канал
	double pole[4][5][2];  //помеха,канал,амплитуда и фаза
						   //double diag[4][5][2];  //помеха,канал,амплитуда и фаза
	double vs[5][2];   //канал,компонента
	double afs[5][2];  //канал,компонента
	double tns[5];    //канал
					  /***COMMENT***/
					  //********************************************************************

					  //for(n=0;n<Nc;n++)CF[n]=SF[n]=0.;
					  //long int uskod[5],uckod[5],urkod[5],uikod[5];
					  //double Ar[(NC-1)*NC];
					  //long double Аr[I*I+I]  -расширенная матрица действительнoй части

					  //int nr=10;//2;     //номер реализации
					  // cout<<"память: " <<coreleft()<<"\n"; getch();
					  //**************************************************************************
					  //**************************************************************************

	FVHT(dt, 60e6, 0, FVST, FVCT, FVpsi);  //ANALCAN  //инициализируем структуры FVST, FVCT, FVpsi
	FVHT(dt, 60.e6, 1, FVST, FVCT, FVpsi); //ANALCANos
	UFVn(FVSN, FVCN); // обнуляем структуры FVSN, FVCN

	UKONn(KN); //обнуляет структуру KN F4D_L.h

	QFO(fmin0, fmax0, 3, fo0, Q);
	KHT(dt, fo0, Q, 1., 1., 0, KT);//инициализируем структуру KT

	QFO(fmin4, fmax4, 1, fo4, Q);
	KHT(dt, fo4, Q, 1., 1.000, 4, KT);

	QFO(fmin5, fmax5, 1, fo5, Q);
	KHT(dt, fo5, Q, 1., 1.000, 5, KT);

	QFO(fmin6, fmax6, 1, fo6, Q);
	KHT(dt, fo6, Q, 1., 1.000, 6, KT);

	QFO(fmin7, fmax7, 1, fo7, Q);
	KHT(dt, fo7, Q, 1., 1.000, 7, KT);

	QFO(fmin8, fmax8, 1, fo8, Q);
	KHT(dt, fo8, Q, 1., 1.000, 8, KT);

	/*
	int nomkon=8;
	cout<< KT.f0[nomkon]<<"  "<<KT.Q[nomkon]<<"  "<<KT.W[nomkon]<<"  "
	<<KT.Ku[nomkon]<<"  "<<KT.g[nomkon]<<"  "<<KT.h0[nomkon]<<"  "
	<<KT.hc[nomkon]<<"  "<<KT.hs[nomkon]<<"  "<<KT.ct[nomkon]<<"  "
	<<KT.st[nomkon]<<"\n";
	getch();
	*/
	QFO(fmin9, fmax9, 1, fo9, Q);
	KHT(dt, fo9, Q, 1., 1.000, 9, KT);


	KALIBROVKA(dt, fmin0, fmax0, 3, kalibr, Nsl, fo0);//????????????????????????????

													  //Nadap=(int)floor(tadap/tzam); //конец адаптации
	Nn = (int)floor(tn / tzam);//начало усреднения сигнала
	Nk = (int)floor(tk / tzam);//конец усреднения сигнала
	cout << "tzam=" << tzam << "  tstep=" << tstep << "\n";
	_getch();

	Ngraf = (int)floor(tgraf*1e9 / tstep + .5);//количество точек на графиках
	int Ngraf1 = (int)(tadap*1e9 / tstep);

	cout << "количество шагов: Ngraf1=" << Ngraf1 << "\n";
	cout << "количество точек: Ngraf=" << Ngraf << "\n";
	_getch();

	double dtgraf = tgraf / (Ngraf - 0);
	cout << "tgraf=" << tgraf << "  dtgraf=" << dtgraf << "\n";

	int Mgraf = (int)(dtgraf / dt + .5);//период вывода на график[нс]
	long int Lgraf = (long int)(tgraf / dt + .5);
	cout << "Mgraf=" << Mgraf << "\n";
	cout << "Lgraf=" << Lgraf << "\n";
	//cout<<"Nw="<<Nw<<"\n";
	cout << "Nadap=" << Nadap << "\n";
	cout << "Nn=" << Nn << "\n";
	cout << "Nk=" << Nk << "\n";
	_getch();

	QFO(fmin1, fmax1, 3, fo1, Q);
	KHT(dt, fo1, Q, 1., 1., 1, KT);

	QFO(fmin2, fmax2, 3, fo2, Q);
	KHT(dt, fo2, Q, 1., 1., 2, KT);

	QFO(fmin3, fmax3, 1, fo3, Q);
	KHT(dt, fo3, Q, 1., 1.078, 3, KT);
	cout << "Q=" << Q << "  fo3=" << fo3 << "\n";
	_getch();

	QFO(fmin4, fmax4, 1, fo4, Q);
	KHT(dt, fo4, Q, 1., 1.000, 4, KT);
	QFO(fmin5, fmax5, 1, fo5, Q);
	KHT(dt, fo5, Q, 1., 1.000, 5, KT);
	QFO(fmin6, fmax6, 1, fo6, Q);
	KHT(dt, fo6, Q, 1., 1.000, 6, KT);

	QFO(fmin7, fmax7, 1, fo7, Q);
	KHT(dt, fo7, Q, 1., 1.000, 7, KT);
	QFO(fmin8, fmax8, 1, fo8, Q);
	KHT(dt, fo8, Q, 1., 1.000, 8, KT);
	QFO(fmin9, fmax9, 1, fo9, Q);
	KHT(dt, fo9, Q, 1., 1.000, 9, KT);
	//-----------------------------------------------------
	UKONn(KN);

	INTCHT(dt, 157e-9, 0, INTCT);
	//INTCHT(dt,230e-9,0,INTCT);
	INTCHT(dt, .06e-6, 1, INTCT);
	//----------------------------------------------
	DIFCHT(dt, 1550e-9, 0, DIFCT);


	ksl = 57439874;//23689732;
	RND(ksl);

	/********************РИСУЕМ ТАБЛИЦУ********************/
	/*************Этот цикл идет до конца main*************/
	int nPRO;
	for (nPRO = 0; nPRO < 1; nPRO++) { //1 BEGIN nPRO
		clrscr();
		window(1, 1, 80, 24);
		gotoxy(65, 1); //gotoxy
		_cprintf("  nPRO=");
		_cprintf("%2d", nPRO);
		gotoxy(1, 1);// gotoxy(1,1)
		_cprintf("Реализации:  Nr=");
		_cprintf("%2d", Nr);
		_cprintf("   nr=");

		gotoxy(30, 1);
		_cprintf("  Шаги:  Nk=");
		_cprintf("%4d", Nk);
		_cprintf("  m=");
		// getch();
		gotoxy(9, 2);
		_cprintf("Загрубление");	//gotoxy
		gotoxy(1, 3);
		_cprintf(" nr     ДО    ПОСЛЕ    Подавл");//gotoxy

		for (n = 0; n < 100; n++) { //600
			X[n] = 0.0;
			Y[n] = 0.0;
		}
		pnormS = 0.;

		//-----------------------------------------------------------
		/*
		// 8 mV
		u_acp1[0][1]=0.008041;  u_acp2[0][1]=0.007935;  u_acp1[1][1]=0.007904;  u_acp2[1][1]=0.007187;
		u_acp1[0][2]=0.007248;  u_acp2[0][2]=0.00766;  u_acp1[1][2]=0.006989;  u_acp2[1][2]=0.006821;
		u_acp1[0][3]=0.007446;  u_acp2[0][3]=0.007507;  u_acp1[1][3]=0.007126;  u_acp2[1][3]=0.006989;
		u_acp1[0][4]=0.008087;  u_acp2[0][4]=0.007584;  u_acp1[1][4]=0.008026;  u_acp2[1][4]=0.00737;
		u_acp1[0][5]=0.007324;  u_acp2[0][5]=0.007782;  u_acp1[1][5]=0.007813;  u_acp2[1][5]=0.00737;
		u_acp1[0][6]=0.007401;  u_acp2[0][6]=0.007584;  u_acp1[1][6]=0.007019;  u_acp2[1][6]=0.007492;
		u_acp1[0][7]=0.007172;  u_acp2[0][7]=0.007614;  u_acp1[1][7]=0.007523;  u_acp2[1][7]=0.00769;
		u_acp1[0][8]=0.007294;  u_acp2[0][8]=0.007477;  u_acp1[1][8]=0.00737;  u_acp2[1][8]=0.008026;
		u_acp1[0][9]=0.007736;  u_acp2[0][9]=0.006989;  u_acp1[1][9]=0.007278;  u_acp2[1][9]=0.007614;
		u_acp1[0][10]=0.007675;  u_acp2[0][10]=0.007126;  u_acp1[1][10]=0.007965;  u_acp2[1][10]=0.007187;
		//	      */
		//-----------------------------------------------------------

		for (nr = 1; nr <= Nr; nr++) { //1 BEGIN Nr
									   //for(nc=0;nc<Nc;nc++) u_acp1[nc][nr]=u_acp2[nc][nr]=0.;

									   //******************************************************************
									   //ВРАЩАЕТСЯ АНТЕННА
			naklon[1] = -45 + 1. / 6.*(nr - 1);   //azimut
			naklon[1] = 0;

			azumant[0][1] = naklon[1];

			for (j = 0; j < Nc; j++)
				for (n = 0; n < 40;n++)  zader[j][n] = 0;

			for (j = 0; j < Nc; j++)
				for (n = 0; n < 25; n++) {
					ucdigzad[j][n] = usdigzad[j][n] = 0;
				}

			//nr=3;
			//   if(nr==11)goto ZZZ17;
			//   azumant[0][1]=-45.+(90./(Nr-1)*(nr-1)); //гл луч

			gotoxy(25, 1);
			_cprintf("%2d", nr); //gotoxy

			gotoxy(2, 3 + nr);
			_cprintf("%2d", nr);//gotoxy

			for (n = 0; n < Nc; n++) Ksl[n] = 1. + .0*RND(ksl);


			UKONn(KN);
			UFVn(FVSN, FVCN);
			UINTCn(INTCN);
			UDIFCn(DIFCN);

			for (n = 0; n < Nc - 1; n++) uwc[n] = uws[n] = Wds[n] = Wdc[n] = .0;
			/* */
			for (n = 0; n < Nc - 1; n++) uints[n] = uintc[n] = .0;

			for (j = 0; j < Nc; j++)
				for (n = 0; n < 40; n++) zader[j][n] = .0;

			for (n = 0; n < Nc;n++) pc[n] = pcc[n] = pcr[n] = pci[n] = .0;
			Pnorm = 0.;

			REALIZ(nr, Np, Nc, v, af, vs, afs, tns);
			ANTCAN(naklon, geom, azumant, pomeha, Np, Nc, zad, pole);

			gotoxy(50, 6);	//gotoxy

			for (np = 0; np < Np; np++) {
				nc = 0;
				_cprintf("  %-3.2f ", 20.*log10(pole[np][nc][0]) - pomeha[np][2]);
			}

			for (j = 0; j < Nc - 1; j++) { Wc[j] = Ws[j] = 0.; }

			L = -1; //mm=0;
			pnorm1 = pn1 = pn2 = pk1 = 0.;
			//goto ZZ456;


			for (m = 0; m < Nk; m++) { //3 BEGIN Nk

									   //takt=66*dt		  
									   //      cout<<"nr="<<nr<<"  m="<<m<<"\n";//getch();
				gotoxy(50, 1); //gotoxy

				_cprintf("%4d\r\n", m);
				//	if(m==142) getch();
				pnorm = pn = pk = 0.;


				//gotoxy(50, 2);//ОТЛАДКА
				//_cprintf("%4d\r\n", Nsl);//ОТЛАДКА

				for (i = 0; i < 2 * Nsl; i++) { //4 BEGIN Nsl
					//gotoxy(50, 2);//ОТЛАДКА
					//_cprintf("i = %d\n", i);//ОТЛАДКА

					L = L + 1;
					t = dt*L;
					//if(m==142&&i==13){cout<<"i="<<i<<"  gdj\n";getch();}
					//if(m==4){cout<<"\nm="<<m<<"  i="<<i<<"\n";getch();}
					//for(j=0;j<Nc;j++)cout<<j<<"  "<<vs[j][0]<<"  "<<vs[j][1]<<"\n";getch();
					//cout<<" Q0 ";

					if (i == 0) {}//


					SIGNAL(dt, fo0, L, Nsl, kalibr, Nc, Np, v, af, tns, vs, afs, zad, pole, usignal);

					//if(L%1000==0)cout<<L<<"  "<<usignal[1]<<"  ";
					//DIVIDER(Nc,0.,usignal);
					//if(L%1000==0){cout<<usignal[1]<<"\n"; getch();}
					//cout<<" Q1 ";
					ZADERG(0, 24, L, usignal[0], usignal[0]);//usignal[0]=uzad;//17нс на AФЧХ REGUS//и наклон полотна
															 //cout<<"m="<<m<<"  i="<<i<<"  usignal0="<<usignal[0]
															 //    <<"  usignal1="<<usignal[1]<<"\n";getch();
					SDVIG(Nc, usignal, usignals, usignalc);

					if (pomeha[2][1]<-10. || pomeha[2][1]>10.) uws[4] = uwc[4] = 0.;

					SUMMATOR1(Nc, usignal, uws, uwc, usigm);
					usigm2 = usigm;

					gotoxy(50, 4);

					//cout<<" Q2 ";


					OKFIL(0, usigm, usigm2);
					//cout<<usigm<<"  "<<usigm2<<"\n";getch();
					//if(nr==2){cout<<"   "<<usignal[0]<<"\n";getch();}
					usigm3 = usignal[0];
					OKFIL(1, usignal[0], usigm3);

					if (t < tadap) {	 // 5 BEGIN t < tadap
										 //usignal[1]=1;
						ZADERG(1, 0, L, usignal[1], usignal[1]);
						ZADERG(2, 0, L, usignal[2], usignal[2]);
						ZADERG(3, 0, L, usignal[3], usignal[3]);
						ZADERG(4, 0, L, usignal[4], usignal[4]);
						ANALCAN(Nc, usignal, us);
						ANALCANos(usigm, usigm1);
						//-----------------------------------------------------
						//
						//cout<<" Q3 ";

						ntakt = L%tstep;
						nstep = (int)floor((float)L / tstep);

						if (ntakt == 21) {
							for (nc = 1; nc < Nc; nc++) {
								ucdig[nc] = us[nc] + .0;
							}//nc
							ucdig[0] = usigm1 + .0;
						}//if(ntakt==21)

						//gotoxy(75, 1);//ОТЛАДКА
						//_cprintf("ntakt = %4d\n", ntakt);//ОТЛАДКА
						//if(ntakt==43)
						if (ntakt == 33) {//33 6 BEGIN 	ntakt==33

							nstep = (int)floor((float)L / tstep);

							for (nc = 1; nc < Nc; nc++) {
								usdig[nc] = us[nc] + .0;
							}//nc

							usdig[0] = usigm1 + .0;

							for (nc = 0; nc < Nc; nc++) {

								q = floor(128.*fabs(ucdig[nc]) + .5);

								if (q > 127) q = 127;

								q = q*SIGN(ucdig[nc]);
								ucdig[nc] = q / 128.;
								q = floor(128.*fabs(usdig[nc]) + .5);

								if (q > 127) q = 127;

								q = q*SIGN(usdig[nc]);
								usdig[nc] = q / 128.;

								//gotoxy(100, 2);//ОТЛАДКА
								//_cprintf("nc = %4d\n", nc);//ОТЛАДКА
								//--------------------------------------------
								/*
								//  Определение постоянной составляющей U АЦП
								u_acp1[nc][nr]=u_acp1[nc][nr]+usdig[nc]/512.;
								u_acp2[nc][nr]=u_acp2[nc][nr]+ucdig[nc]/512.;
								*/
								//--------------------------------------------
								/*
								//  Компенсация постоянной составляющей U АЦП
								//cout<<L<<"aaaaaaaaa"<<usdig[nc]<<"sssssssssss";
								if(nc<2)m_f=nc;
								else m_f=1;
								usdig[m_f]=usdig[m_f]-u_acp1[m_f][nr];
								//cout<<usdig[nc]<<"\n"; getch();
								ucdig[m_f]=ucdig[m_f]-u_acp2[m_f][nr];
								*/
								//--------------------------------------------

								//usdig[nc]=floor(usdig[nc]*128.)/128.;
								//cout<<nstep<<"  "<<  usdig[nc]<<"  "<<  ucdig[nc]<<"\n";

							} //for(nc=0; nc<Nc; nc++)

							  //getch();

							  //if(nstep<480)
							  //cout<<"Nstep="<<Nstep<<"  nstep="<<nstep<<"  ntakt="<<ntakt<<"\n";getch();
							  //       if(ntakt==33&&nstep<Nstep)
							  //{PROCESSORdig(Nc, Np,nstep,usdig, ucdig, uints, uintc,Wds,Wdc,Pnorm);
							PROCESSORdig1(Nc, Np, nstep, usdig, ucdig, uints, uintc, Wds, Wdc, Pnorm);
							//Dpzu<<nstep<<"  "<<(pzu*1.)<<"\n"; // cout<<"processor\n";getch();

							// PROCESSORdig(Nc, Np,nstep,usdig, ucdig, uints, uintc,Wds,Wdc);
							//PROCESSOR(Nc, nstep,usdig, ucdig, uints, uintc,    nPRO);
							//PROCESSOR_12(Nc, nstep,usdig, ucdig, uints, uintc);
							//PROCESSOR_14(Nc, nstep,usdig, ucdig, uints, uintc);
							// PROCESSOR_4(Nc,nstep,pnormS,us,uc,uints, uintc,pc,pcr,pci,kor);
							//PROCESSOR_STRL(Nc, nstep, usdig, ucdig,uints, uintc, Wds, Wdc);
							//--------------------------------------------------------------------
							//                             Задержка конвейера
							zadacp = 0;
							if (zadacp > 0) {
								jjj = nstep%zadacp;
								for (nc = 0; nc < Nc; nc++) {

									a = usdigzad[nc][jjj];
									usdigzad[nc][jjj] = uints[nc];
									uints_zad[nc] = a;

									a = ucdigzad[nc][jjj];
									ucdigzad[nc][jjj] = uintc[nc];
									uintc_zad[nc] = a;
								}//for nc
							}
							else {

								for (nc = 0; nc < Nc; nc++) {
									uints_zad[nc] = uints[nc];
									uintc_zad[nc] = uintc[nc];
								} //for nc
							} //else

							  //--------------------------------------------------------------------
							  //}if(nstep<480)
							  //cout<<"nstep="<<nstep<<"  "<<uintc[0]<<"  ";

							DTkk << nstep << "  ";
							for (n = 1; n < Nc; n++) DTkk << kor[n] << "  ";
							DTkk << "\n";
							//cout<<"  step="<<floor(L/210.)<<"  ntakt="<<ntakt
							//    <<"  "<<us[0]<<"  "<<uc[0]<<"  "<<"  "<<us[1]<<"  "<<uc[1]
							//     <<"  "<<uints[1]<<"  "<<uintc[1]<<"\n";getch();

						}//6 END if(ntakt==33)
					} // 5 END if(t<tadap)

					  //uintc[0]=-.284;  uints[0]=.500; //ручной режим

					REGUS(L, Nc, usignals, usignalc, uints_zad, uintc_zad, uws, uwc, Wc, Ws);
					//if(m>700)
					//{cout<<"L="<<L<<"  m="<<m<<"  "<<uints[0]<<"  "<<uintc[0]
					//    <<"  Wc[0]="<<Wc[0]<<"  Ws[0]="<<Ws[0]<<"\n"; getch();
					// }

					//*********************************************************

					//                     загрубление
					pn = pn + usigm3*usigm3;
					pk = pk + usigm2*usigm2;
					pnorm = pnorm + ashum*ashum*.5;

					//pk1=pk1+usigm1*usigm1;
					if (m > Nn)//if(m>Nust+Ngraf)
					{
						pn1 = pn1 + usigm3*usigm3;
						//    pn2=pn2+us[0]*us[0];
						pn2 = pn2 + usignal[1] * usignal[1];
						//    pn2=pn2+usigm1*usigm1;
						pk1 = pk1 + usigm2*usigm2;
						pnorm1 = pnorm1 + ashum*ashum*.5;
					}//if(m > Nn)

					 //*****************************************************************
					if (L < Lgraf && L % Mgraf == 0) { //5 BEGIN if((int)(t*1e9)%(int)(tzam*1e9)==0&&t<=tgraf)
						podav = DB(pn / pnorm);
						podav1 = DB(pk / pnorm);
						//X[L / Mgraf] = X[L / Mgraf] + (float)(pn / pnorm);
						//Y[L / Mgraf] = Y[L / Mgraf] + (float)(pk / pnorm);

						if (m == 0) {
							DTw << "\n   " << Ngraf << "  " << dtgraf << "  " << Nc - 1 << "\n";
							DTw_dig << Ngraf << "  " << dtgraf << "  " << Nc - 1 << "\n";
							DTz << Ngraf << "  " << dtgraf << "  " << Nc - 1 << "\n";
							DTkk << Ngraf1 << "  " << Nc - 1 << "\n";

							if (nr == 1)  DTzs << Ngraf << "  " << dtgraf << "  " << Nc - 1 << "\n";

						} //if(m==0)

						  //DTw<<m<<"  ";
						DTw << L / Mgraf << "  ";
						DTw_dig << L / Mgraf << "  ";

						for (j = 0; j < Nc - 1; j++) { DTw << (float)Wc[j] << "  " << (float)Ws[j] << "  "; }

						DTw << "\n";

						for (j = 0; j < Nc - 1; j++) { DTw_dig << (float)uintc[j] << "  " << (float)uints[j] << "  "; }

						DTw_dig << "\n";

						DTz << L / Mgraf << "  " << L*dt*1e6 << "  " << podav << "  " << podav1;
						DTz << "\n";
					}// 5 END if(L < Lgraf && L%Mgraf == 0)

					 //*********************************************************
					 //*********************************************************

				}//4 END Nsl

				 //*********************************************************
				podav = DB(pn / pnorm);
				podav1 = DB(pk / pnorm);
				//*********************************************************
			}//3 END Nk
			 //*********************************************************
		ZZ456:

			podav2 = DB(pn1 / pnorm1);
			podav21 = DB(pn2 / pnorm1);
			podav3 = DB(pk1 / pnorm1);
			DTz << "  " << podav2 << "  " << podav3 << "\n";
			DTze << nr << "  " << podav2 << "  " << podav3 << "  " << podav2 - podav3 << "\n";
			DTzsa << nr << "  " << podav2 << "  " << podav3 << "  " << podav2 - podav3 << "\n";
			//      DTp<<nr<<"  "<<pomeha[0][1]<<"  "<<podav2<<"  "<<podav3<<"  "
			//	   <<podav2-podav3<<"\n";
			DTp << nr << "  " << naklon[1] << "  " << podav2 << "  " << podav3 << "  " << podav2 - podav3 << "\n";
			DTwe << nr << " ";
			DTwsa << nr << " ";

			for (j = 0;j < Nc - 1;j++) DTwe << Wc[j] << "  " << Ws[j] << "  ";DTwe << "\n";
			for (j = 0;j < Nc - 1;j++) DTwsa << Wc[j] << "  " << Ws[j] << "  ";DTwsa << "\n";
			zagr[nr - 1] = (float)podav3;
			zagr1[nr - 1] = (float)podav2;
			//*********************************************************

			//******************************************************************8
			/*
			//ВРАЩЕНИЕ ПОЛОТНА

			ww=sqrt(Wc[0]*Wc[0]+Ws[0]*Ws[0]);
			DTdiag<<"\n"<<nr<<"  ";
			for(pp=0;pp<6;pp++)
			{ naklon[1]=naklon[1]+.1*pp;
			azumant[0][1]=azumant[0][1]+.1*pp;
			azumant[1][1]=azumant[1][1]+.1*pp;
			ANTCAN_1(naklon,geom,azumant,pomeha,Np,Nc,zad,diag);
			naklon[1]=naklon[1]-.1*pp;
			azumant[0][1]=azumant[0][1]-.1*pp;
			azumant[1][1]=azumant[1][1]-.1*pp;
			re0=diag[0][0][0]*cos(diag[0][0][1]);im0=diag[0][0][0]*sin(diag[0][0][1]);
			re1=diag[0][1][0]*cos(diag[0][1][1]);im1=diag[0][1][0]*sin(diag[0][1][1]);
			psiww=diag[0][0][1]-diag[0][1][1];
			if(pp==0){wwc=ww*cos(psiww);wws=ww*sin(psiww);}
			//res=re0+(Wc[0]*re1-Ws[0]*im1); ims=im0+(Wc[0]*im1+Ws[0]*re1);
			diag3c=re0-re1*wwc;diag3s=im0-im1*wws;
			diag3=(float)(10*log10(diag3c*diag3c+diag3s*diag3s));
			diag1=(float)(10*log10(re0*re0+im0*im0));
			diag2=(float)(10*log10(re1*re1+im1*im1));
			if(pp==0)DTdiag<<diag1<<"  ";
			DTdiag<<diag3<<"  ";
			//  DTdiag<<diag1<<"  "<<diag2<<"  "<<diag3<<"  "<<podav2<<"  "<<podav3<<"\n";
			}
			*/
			//******************************************************************8
			/*
			cout<<"ЗАГРУБЛЕНИЕ: без адаптации-  "<<podav2<<" db\n"<<
			"             c адаптацией-   "<<podav3<<" db\n"<<
			"ПОДАВЛЕНИЕ:  "<<podav2-podav3<<" db\n";
			cout<<"dWmax="<<dW<<"\n";
			*/
			gotoxy(6, 3 + nr);  //gotoxy
			_cprintf("podav2 %5.3f ", podav2);
			_cprintf(" %5.3f ", podav3);
			_cprintf("  %5.3f ", podav2 - podav3);
			_cprintf("  %5.3f ", podav21);
			//getch();getch();
		ZZZZ545:
			if (nr == Nr) {
				for (n = 0; n < 600; n++) {
					X[n] = (float)(10.*log10(X[n] / Nr + 1e-12));
					Y[n] = (float)(10.*log10(Y[n] / Nr + 1e-12));

					DTzs << n << "  " << n*dt*1e6 << "  " << X[n] << "  " << Y[n] << "\n";
				}//n
			}//if(nr==Nr)
			 //getch();
			 //-------------------------------------------------------
			 /*
			 DTukor <<nr<<"  ";
			 for(nc=0;nc<Nc;nc++)    DTukor<<"u_acp1["<<nc<<"]["<<nr<<"]="
			 <<u_acp1[nc][nr]
			 <<"  u_acp2["<<nc<<"]["<<nr<<"]="
			 <<u_acp2[nc][nr]
			 <<"  ";
			 DTukor <<"\n";
			 */
			 //-------------------------------------------------------

		} //2 END Nr 

		DTukor << "--------------------------------\n\n";
	ZZZ17:
	ZZZZ546:

		// for(n=0;n<Nr;n++)cout<<"nr="<<n<<"  загрубление="<<zagr[n]<<"db\n";
		textattr(4 + (7 << 4));
		gotoxy(1, 3 + Nr + 1); //gotoxy
		_cprintf("средн");

		podav2 = 0.;
		for (n = 0;n < Nr;n++)podav2 = podav2 + pow(10., zagr1[n] / 10.);
		podav2 = DB(podav2 / Nr);
		gotoxy(7, 3 + Nr + 1);//gotoxy
		_cprintf("%5.3f  ", podav2);
		//cout<<"среднее загрубление без адаптации="<<podav2<<"db\n";

		podav3 = 0.;
		for (n = 0; n < Nr; n++)   podav3 = podav3 + pow(10., zagr[n] / 10.);
		podav3 = DB(podav3 / Nr);
		gotoxy(15, 3 + Nr + 1);//gotoxy
		_cprintf("%5.3f", podav3);
		gotoxy(24, 3 + Nr + 1); //
		_cprintf("%5.3f \r\n", podav2 - podav3);
		//cout<<"среднее загрубление с адаптацией="<<podav3<<"db\n";

		DTze << podav2 << "  " << podav3 << "  " << podav2 - podav3 << "\n";
		DTzsa << podav2 << "  " << podav3 << "  " << podav2 - podav3 << "\n";
		// sound(1000);delay(500);sound(10);
		//getch();getch();

	}//1 END nPRO 

	_getch();
	_getch();
NM:;
	//1472 
	return 0;
}

//*************************************************

void QFO(float f1, float f2, int n, float& fo, float& Q) {
	//f1 - Fнижн; f2- Fверх; n- кoл-вo кoнтурoв
	fo = (float)sqrt(f1*f2);
	//Q=(float)pow(.5,1./n);
	Q = sqrt(pow(2., 1. / n) - 1.);
	//if((1-Q)/Q<=0.){cout<<"под корнем QFO меньше 0\n";getch();getch();}
	Q = (float)(fo / (f2 - f1)*Q);
}  //QFO

   //**************************************************************

void FVHT(float dt, float f0, int nomfvtip, FVTIP& FVST, FVTIP& FVCT, FVp& FVpsi) {
	float as, bs, ac, bc, s, c, h;
	FVpsi.f0[nomfvtip] = f0;  //средняя частота диапазона 60e6(Гц)
	FVpsi.h[nomfvtip] = h = FVpsi.f0[nomfvtip] / 2550.;
	FVpsi.as[nomfvtip] = as = 770.*h;
	FVpsi.bs[nomfvtip] = bs = 8330.*h;
	FVpsi.ac[nomfvtip] = ac = 3080.*h;
	FVpsi.bc[nomfvtip] = bc = 33320.*h;
	FVST.a[nomfvtip] = exp(-dt*as);
	FVST.b[nomfvtip] = exp(-dt*bs);
	FVCT.a[nomfvtip] = exp(-dt*ac);
	FVCT.b[nomfvtip] = exp(-dt*bc);
	s = 2.*(as + bs) / ((as - bs)*dt);
	c = 2.*(ac + bc) / ((ac - bc)*dt);
	FVST.ga[nomfvtip] = s*(1. - FVST.a[nomfvtip]) / as;
	FVST.gb[nomfvtip] = s*(1. - FVST.b[nomfvtip]) / bs;
	FVCT.ga[nomfvtip] = c*(1. - FVCT.a[nomfvtip]) / ac;
	FVCT.gb[nomfvtip] = c*(1. - FVCT.b[nomfvtip]) / bc;
	FVST.k[nomfvtip] = FVCT.k[nomfvtip] = .5;
	//  cout<<"Параметры фазосдвигателя:\nas="<<as<<"  bs="<<bs<<"  ac="<<ac<<
	//	"  bc="<<bc<<"\n"; getch();getch();
} //FVHT

  //***************************************************************************
void UFVn(FVNOM& FVSN, FVNOM& FVCN) {
	int n;
	for (n = 0; n < Nfv; n++)
		FVSN.u0[n] = FVSN.ua[n] = FVSN.ub[n] = FVSN.du[n] = FVSN.un[n] =
		FVCN.u0[n] = FVCN.ua[n] = FVCN.ub[n] = 0.;
}  //	UFVn

   //*************************************************
void KHT(float dt, float f0, float Q, float W, float Ku, int nomkontip, KONTIP& KT) {

	float c, c1, del, hi, ex, dd;
	KT.Q[nomkontip] = Q;//=6.;//12.; //добротность
	KT.W[nomkontip] = W;//=1.; //волновое сопротивление контура
	KT.f0[nomkontip] = f0;//=6e7;//*1e-9;  f0=1/(2*PI*sqrt(L*C)) (Гц)
	KT.Ku[nomkontip] = Ku;//=.5;  //усиление по ннапряжению на резонанс. частоте
	KT.g[nomkontip] = Ku / (W*Q);
	c = 2.*M_PI*f0; //круговая частота омега0
					//cout<<"tip "<<nomkontip<<"  fo="<<f0<<"\n";getch();
	c1 = c*(float)sqrt(1. - .25 / (Q*Q));// омега рез.
	del = c*.5 / Q*1.;  //затухание дельта
	hi = atan(1. / (Q - .5)) + atan(1. / (2 * Q - .5 / Q));
	ex = exp(-del*dt);
	dd = W*(float)sqrt(1. - 1 / Q + 1.25 / (Q*Q)) / (c*dt);
	KT.hc[nomkontip] = dd*(cos(hi) - ex*cos(c1*dt - hi));
	KT.hs[nomkontip] = dd*(sin(-hi) - ex*sin(c1*dt - hi));
	KT.h0[nomkontip] = W / Q;
	KT.ct[nomkontip] = ex*cos(c1*dt);KT.st[nomkontip] = ex*sin(c1*dt);
} //KHT

  //******************************************************************
void DIFCHT(float dt, float tau, int tip, DIFCTIP& DIFCT) {
	DIFCT.tau[tip] = tau / dt; //tau=3e-7[сек] -постоянная времени диффер.цепочки
	DIFCT.e[tip] = exp(-dt / tau);
}//FNCHT

 //***************************************************************************
void UDIFCn(DIFCNOM& DIFCN) {
	int n;
	for (n = 0; n < Ndifc; n++)
		DIFCN.du[n] = DIFCN.u0[n] = DIFCN.u1[n] = DIFCN.un[n] = 0.;
} //UFNCn


  //********************************************************
void ANTCAN(float naklon[2], float geom[5][2], float azumant[5][3], float pomeha[4][3], int Np, int Nc, double zad[4][5], double pole[4][5][2])
{
	int np, nc, m, n;
	double re, im, st, sf, cf, lmin;
	float tefi[4][2]; //[номер помехи][коорд помехи]
	float tefi0[5][2]; //[номер антенны][коорд гл луча]
	lmin = 100;
	for (np = 0; np < Np; np++) {
		TEFI(naklon[0], naklon[1], pomeha[np][0], pomeha[np][1], tefi[np][0], tefi[np][1]);
		//cout<<"np="<<np<<"  tet="<<tefi[np][0]<<"  fi="<<tefi[np][1]<<"\n";
		//getch();
		for (nc = 0;nc < Nc;nc++) {
			TEFI(naklon[0], naklon[1], azumant[nc][0], azumant[nc][1], tefi0[nc][0], tefi0[nc][1]);
			// if(nc==0)
			D400(tefi0[nc][0] * GR, tefi0[nc][1] * GR, tefi[np][0] * GR, tefi[np][1] * GR,
				normir, azumant[nc][2], amp, re, im);
			//cout<<"re="<<re<<"  im="<<im<<"\n";
			if (re*re + im*im == 0.) {
				cout << "под корнем ANTCAN меньше 0\n";
				_getch();
				_getch();
			}
			pole[np][nc][0] = sqrt(re*re + im*im)*
				pow(10., (pomeha[np][2] + G[(int)azumant[nc][2]]) / 20.);//амп
																		 //cout<<"np="<<np<<"  nc="<<nc<<"  pole="<<pole[np][nc][0]<<"\n";getch();
																		 //cout<<"np="<<np<<"  nc="<<nc<<" azumant[nc][2]="<<azumant[nc][2]
																		 //    <<"  "<<G[(int)azumant[nc][2]]<<"\n";getch();
			pole[np][nc][1] = ARCTG(re, im);//фаза
			st = sin(tefi[np][0] * GR);
			sf = sin(tefi[np][1] * GR);
			cf = cos(tefi[np][1] * GR);
			zad[np][nc] = -st*(geom[nc][0] * cf + geom[nc][1] * sf);//в длинах волн
			if (zad[np][nc] < lmin) lmin = zad[np][nc];
			//cout<<"zad="<<zad[np][nc]<<"  lmin="<<lmin<<"\n";getch();
		}//nc
	}//np
	 //cout<<"lmin="<<lmin<<"\n";getch();
	for (np = 0; np < Np; np++)
		for (nc = 0; nc < Nc; nc++) {
			zad[np][nc] = zad[np][nc] - lmin;
			pole[np][nc][1] = pole[np][nc][1] - 2 * M_PI*zad[np][nc];//фаза
																	 //cout<<"zad="<<zad[np][nc]<<"  lmin="<<lmin<<"\n";//getch();
			zad[np][nc] = zad[np][nc] * lam / 3e11;//3e11- скор света[mм/сек]
		}//nc
		 //cout<<"ANTCAN:   lam="<<lam<<"\n";getch();
}//ANTCAN

 //********************************************************
void SIGNAL(double dt, double fo, long int n, int Nsl, double kalibr, int Nc, int Np, double v[4][5][2], double af[4][5][2],
	double tns[5], double vs[5][2], double afs[5][2],
	double zad[4][5], double pole[4][5][2],
	double usignal[5])
{
	int nc, np, m, q;       //double uv_f,u_f;
	int nomkontip = 0;
	//cout<<"shum  ";
	for (nc = 0;nc < Nc;nc++) {
		//cout<<afs[nc][0]<<"  "<<afs[nc][1]<<"  "<<usignal[nc]<<"  ";
		//cout<<"signal_0"<<"\n";
		SHUM(dt, fo, n, Nsl, tns[nc], vs[nc], afs[nc], u_f);
		usignal[nc] = u_f;
		//cout<<nc<<vs[nc][0]<<"  "<<vs[nc][1]<<"  "<<afs[nc][0]<<"  "<<afs[nc][1]<<"  "<<usignal[nc]<<"  ";
		for (np = 0; np < Np; np++) {
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

 //********************************************************
void UPOMCAN(double dt, double fo, long int n, int N, double zad, double v[2], double af[2], double pole[2], double& u) {//         v-  2 аргумента для RALLEY
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

 //**************************************************************
void SDVIG(int Nc, double usignal[], double usignals[], double usignalc[]) {
	int n_f;
	for (n_f = 1; n_f < Nc; n_f++)
	{
		UFV(n_f, 0, FVST, FVCT, FVSN, FVCN, usignal[n_f], usignals[n_f], usignalc[n_f]);
		usignals[n_f] = 2 * usignals[n_f];
		usignalc[n_f] = 2 * usignalc[n_f];
	}//n
}//SDVIG

 //***************************************************************************
void UFV(int n, int tip, FVTIP&FVST, FVTIP&FVCT, FVNOM&FVSN, FVNOM&FVCN,
	double uv, double&us, double&uc) {

	FVSN.du[n] = uv - FVSN.un[n];
	FVSN.u0[n] = FVSN.u0[n] + FVSN.du[n];
	FVSN.ua[n] = FVSN.ua[n] * FVST.a[tip] + FVSN.du[n] * FVST.ga[tip];
	FVSN.ub[n] = FVSN.ub[n] * FVST.b[tip] - FVSN.du[n] * FVST.gb[tip];
	us = (FVSN.u0[n] + FVSN.ua[n] + FVSN.ub[n])*FVST.k[tip];
	FVCN.u0[n] = FVCN.u0[n] + FVSN.du[n];
	FVCN.ua[n] = FVCN.ua[n] * FVCT.a[tip] + FVSN.du[n] * FVCT.ga[tip];
	FVCN.ub[n] = FVCN.ub[n] * FVCT.b[tip] - FVSN.du[n] * FVCT.gb[tip];
	uc = (FVCN.u0[n] + FVCN.ua[n] + FVCN.ub[n])*FVCT.k[tip];
	FVSN.un[n] = uv;
} //UFV

  //*******************************************************
void ZADERG(int nom, int nzad, long int n, double u1, double& u2) {
	if (nzad > 0) {
		u2 = zader[nom][n % (nzad)];
		zader[nom][n % (nzad)] = u1;
	}
	else u2 = u1;
}// ZADERG

 //********************************************************
void SHUM(double dt, double fo, long int n, int Nsl, double tns, double vs[2], double afs[2], double& u) {
	//u- нормир шум
	int m;

	m = (int)floor(tns*Nsl);
	//cout<<"m+1="<<m+1<<"  n%Nsl="<<n%Nsl<<"\n";
	if (n%Nsl == (m + 1) % Nsl) {
		RALLEY(vs, afs);
		if (n < Nsl)
			RALLEY(vs, afs);
		//cout<<"\ntjrop\n";
	}

	u = afs[0] * sin(2 * M_PI*fo*(n - m)*dt + afs[1]);

}//SHUM

 //********************************************************
void REGUS(int L, int Nc, double usignals[], double usignalc[], double uints[], double uintc[], double uws[], double uwc[], double Wc[], double Ws[]) {
	int m_f, n_f, kontip, intctip; //double u1,u2;
	float a = 0.;   //минимальное усиление
	double b, c, e, f;
	//33<=nkon<41
	kontip = 3;
	intctip = 0;
	for (n_f = 0; n_f < Nc - 1; n_f++) {
		m_f = 6 * n_f;
		UINTC(m_f + 0, intctip, INTCT, INTCN, uints[n_f], u_f);
		UINTC(m_f + 1, intctip, INTCT, INTCN, u_f, uv_f);
		UINTC(m_f + 2, intctip, INTCT, INTCN, uv_f, Ws[n_f]);
		//     Ws[n_f]=Ws[n_f]-.9;
		//Ws[n_f]=Ws[n_f]*.5*1.7/6.25;

		//Ws[n_f]=.588*uints[n_f]*4;

		UINTC(m_f + 3, intctip, INTCT, INTCN, uintc[n_f], u_f);
		UINTC(m_f + 4, intctip, INTCT, INTCN, u_f, uv_f);
		UINTC(m_f + 5, intctip, INTCT, INTCN, uv_f, Wc[n_f]);
		//       Wc[n_f]=Wc[n_f]+.9;
		//Wc[n_f]=Wc[n_f]*.5*1.7/6.25;

		//Wc[n_f]=.588*uintc[n_f]*4;

		//------------------------------------------------------
		/*
		//просачивание управляющего напряжения на сигнальный вход
		UDIFC(2*n_f,0,DIFCT,DIFCN,Ws[n_f],u_f);
		usignals[n_f+1]=usignals[n_f+1]+.0000000316*u_f;
		UDIFC(2*n_f+1,0,DIFCT,DIFCN,Ws[n_f],u_f);
		usignalc[n_f+1]=usignalc[n_f+1]+.0000000316*u_f;
		*/
		//------------------------------------------------------
		//собственный шум REGUS
		b = .0000088*sin(L*M_PI * 2 * (.060 + n_f*.0003))*1.1;
		c = .0000088*cos(L*M_PI * 2 * (.060 + (n_f - 4)*.0003))*1.1;
		//  b=c=0.;  //обнуление собственного шума REGUS
		//cout<<"L="<<L<<"  "<<"b="<<b<<"\n";getch();
		if (fabs(Ws[n_f]) < a)  e = SIGN(Ws[n_f])*a;  else e = Ws[n_f];
		if (fabs(Wc[n_f]) < a)  f = SIGN(Wc[n_f])*a;  else f = Wc[n_f];

		usignals[n_f + 1] = usignals[n_f + 1] * e + b;
		usignalc[n_f + 1] = usignalc[n_f + 1] * f + c;
		//------------------------------------------------------
		//частотная характеристика рег.ус
		//UKON(33+2*n_f,kontip,KT,KN,(usignals[n_f+1])*
		//	       (a*SIGN(Ws[n_f])+(1-a)*Ws[n_f])+b,uws[n_f]);
		//UKON(33+2*n_f+1,kontip,KT,KN,usignalc[n_f+1]*
		//	       (a*SIGN(Wc[n_f])+(1-a)*Wc[n_f])+c,uwc[n_f]);

		UKON(33 + 2 * n_f, kontip, KT, KN, usignals[n_f + 1], uws[n_f]);
		UKON(33 + 2 * n_f + 1, kontip, KT, KN, usignalc[n_f + 1], uwc[n_f]);
		//uws[n_f]=uws[n_f]+12e-6*sin(L*PI*2*(60+.3*n_f));
		//uwc[n_f]=uwc[n_f]+12e-6*sin(L*2*PI*(60-.3*n_f));

		//----------------------------------------------------
		/*
		//"идеальная" частотная хар-ка рег ус
		uwc[n_f]=usignalc[n_f+1]*Wc[n_f];
		uws[n_f]=usignals[n_f+1]*Ws[n_f];
		*/
		//---------------------------------------------------

	}//n_f
}//REGUS


 //********************************************************
void PROCESSORdig1(int Nc, int Np, int nstep, double udigs[], double udigc[], double ws[], double wc[], double Wds[], double Wdc[], double Pnorm) {
	int nc;
	float pnorm[4], mronp[4], ds, dc, dmax, dd, da, db, fi, cfi, sfi, pzu;
	long int q, q1;
	dmax = 1. / 8.;;

	//dmax=1.;
	//pnorm=0.0003; //добавка для выбора постоянной составляющей на входе АЦП
	// pnorm=0;
	for (nc = 1;nc < Nc;nc++) {
		//pnorm=pnorm+(udigs[nc]*udigs[nc]+udigc[nc]*udigc[nc]);
		pnorm[nc] = (udigs[nc] * udigs[nc] + udigc[nc] * udigc[nc]);
	}
	//pnorm=0.075; //нормировка для помехи 47 дБ
	//pnorm=0.075/3.16; //нормировка для помехи 42 дБ
	//pnorm=0.075/3.16/5.; //нормировка для помехи 35 дБ
	//pnorm=.022/5.;
	//---------------------------------------------------------------------
	/*
	if(nstep<128)Pnorm=Pnorm+pnorm;
	if(nstep==128)Pnorm=Pnorm/128;
	if(nstep>=128)pnorm=Pnorm;       //   усреднение в начале
	//                        */
	//------------------------------------------------------------------------

	// pnorm=.28;
	// pnorm=.85;
	//  pzu=.019*exp(-(nstep)/500.)*2.75*20; //1 к  1 п  при помехе 34db с огр
	//  pzu=.019*exp(-(nstep)/500.);    //1 канал  1 помеха  при помехе 54db
	//  pzu=.019*exp(-(nstep)/500.)*2.75; //1 к  1 п  при помехе 54db  с огр
	//  pzu=.019*exp(-(nstep)/500.);    //1 канал  1 помеха  при помехе 42db
	//pzu=.30*exp(-(nstep)/500.);    // 4 канала  3 помехи:-2,2,3
	//pzu=.19*exp(-(nstep)/500.);    // 4 канала  3 помехи:2,3,4
	//pzu=.12*exp(-(nstep)/500.);    // 4 канала  2 помехи:-2,2
	//pzu=.019*exp(-(nstep)/500.);    // 4 канала  1 помехи:-2,2

	//cout<<"nstep="<<nstep<<"\n";getch();

	ds = exp(-(nstep) / 500.);
	ds = floor(1024.*ds + 1.) / 1024.;
	if (Np == 3)  pzu = .09*ds * 3;
	if (Np == 2)  pzu = .03*ds;// *2.4;//  /2.;
	if (Np == 1)  pzu = .022*ds;
	//	      if(Np==2)pzu=.05*ds;
	//pzu=.5;
	//pzu=.3*exp(-(nstep)/500.);    //
	// pzu=.5*exp(-(nstep)/10000.)*1.;    //
	//pzu=.017*exp((nstep)/900.)*1.;    //
	//	pzu=.20;
	// pzu=.022*1*exp(-(nstep)/500.);
	//pzu=1.7*exp(-(nstep)/90.);       //1000ns

	pzu = floor(pzu*128. + .5) / 128.;

	for (nc = 1;nc < Nc;nc++) {

		if (pnorm[nc] == 0)  mronp[nc] = 0.;  else mronp[nc] = pzu / pnorm[nc];  //7.727 db
	}
	//----------------------------------------------------------------------
	/*
	clrscr();
	cout<<"nstep="<<nstep<<"\n";
	for(nc=0;nc<Nc;nc++)
	{cout<<"       us="<<udigs[nc]*128<<"  uc="<<udigc[nc]*128<<"\n";}
	cout<<"               pnorm="<<pnorm*16384<<"  mronp="<<mronp<<"\n";getch();
	*/
	//----------------------------------------------------------------------
	for (nc = 1;nc < Nc;nc++)
	{
		ds = (udigs[nc] * udigc[0] - udigc[nc] * udigs[0])*mronp[nc];
		dc = (udigc[nc] * udigc[0] + udigs[nc] * udigs[0])*mronp[nc];

		//-------------------------------------------
		//	   //dc=dc*lam; //q=floor(16384.*dc*4.);
		//		    q=(long int)(SIGN(dc)*floor(16384.*fabs(dc)*4.));
		//q=(long int)(floor(fabs(dc)*16384.*4.));
		//		      q=q%(2048*4);dc=q/16384./4.;
		//	   //ds=ds*lam; //q=floor(16384.*ds*4.);
		//		      q=SIGN(ds)*floor(16384.*fabs(ds)*4.);
		//		      q=q%(2048*4);ds=q/16384./4.;
		//-------------------------------------------
		fi = 3.1415 / 180.*300.;  //с задержкой 23ns в компенсационных каналах
		fi = 3.1415 / 180.*20.;  //с задержкой 4ns в компенсационных каналах
		fi = 3.1415 / 180.*75.;  //без задержки в компенсационных каналах
								 //fi=3.1415/180.*115.;
		cfi = cos(fi);      sfi = sin(fi);
		dd = dc;
		dc = dd*cfi + ds*sfi;
		ds = ds*cfi - dd*sfi;
		dd = fabs(dc);
		if (dd > dmax) dc = SIGN(dc)*dmax;
		dd = fabs(ds);
		if (dd > dmax) ds = SIGN(ds)*dmax;

		ws[nc - 1] = ws[nc - 1] + ds;
		//if(nc==1){clrscr();cout<<nstep<<"  "<<ws[nc-1];}
		dd = fabs(ws[nc - 1]);
		if (dd > 1.) ws[nc - 1] = SIGN(ws[nc - 1]);
		else ws[nc - 1] = SIGN(ws[nc - 1])*floor(dd*2048. + .5) / 2048.;

		wc[nc - 1] = wc[nc - 1] + dc;
		//if(nc==1){cout<<"  "<<ws[nc-1]<<"    "<<wc[nc-1];}
		dd = fabs(wc[nc - 1]);
		if (dd > 1.) wc[nc - 1] = SIGN(wc[nc - 1]);
		else wc[nc - 1] = SIGN(wc[nc - 1])*floor(dd*2048. + .5) / 2048.;
	} //nc
	  //--------------------------------------------------------------
	  /*
	  //    		усреднение в конце
	  if(nstep>=Nstep-256&&nstep<Nstep)
	  {  for(nc=1;nc<Nc;nc++)
	  {Wdc[nc-1]=Wdc[nc-1]+wc[nc-1];
	  Wds[nc-1]=Wds[nc-1]+ws[nc-1];
	  }
	  //cout<<"NSTEP="<<nstep<<"\n";  getch();
	  }//if(nstep>480-32&&nstep<480)
	  if(nstep==Nstep-1)
	  {//cout<<"NSTEP  ="<<nstep<<"\n";  getch();

	  for(nc=1;nc<Nc;nc++)
	  {
	  //cout<<"wc["<<nc-1<<"]="<<wc[nc-1]<<"   \n"
	  //    <<"ws["<<nc-1<<"]="<<ws[nc-1]<<"   \n";
	  wc[nc-1]=Wdc[nc-1]/256.;
	  ws[nc-1]=Wds[nc-1]/256.;
	  //cout<<"Wdc["<<nc-1<<"]="<<Wdc[nc-1]<<"   \n"
	  //    <<"wc["<<nc-1<<"]="<<wc[nc-1]<<"   \n"
	  //    <<"Wds["<<nc-1<<"]="<<Wds[nc-1]<<"   \n"
	  //    <<"ws["<<nc-1<<"]="<<ws[nc-1]<<"   \n";         getch();getch();
	  }
	  }//if(nstep>480-32&&nstep<480)
	  //			  */
	  //----------------------------------------------------------------


	  //     ws[nc-1]=Wds[nc-1];
	  //     wc[nc-1]=Wdc[nc-1];
	  //if(nc==1){cout<<"  "<<wc[nc-1]<<"\n";getch();}
	  //cout<<" nc="<<nc<<"  ws="<<ws[nc]<<"  wc="<<wc[nc]
	  //      <<"  ds="<<ds<<"  wc="<<dc<<"\n";   getch();
}//PROCESSORdig1

 //**************************************************************
void ANALCAN(int Nc, double signal[], double udig[]) {
	//усиление основного канала-200 , компенсационного-300
	//long int rty=7;
	// double uv_f,u_f;

	int m_f, n_f;
	int nomkontip = 1;          //  Nc=2;
	for (n_f = 1; n_f < Nc; n_f++) {

		uv_f = (signal[n_f]);
		for (m_f = 0; m_f < 3; m_f++) {
			UKON(15 + n_f * 3 + m_f, nomkontip, KT, KN, uv_f, u_f);
			uv_f = u_f;
		}//m
		uv_f = uv_f*.625*355.;
		udig[n_f] = uv_f;

	}//n
}//ANALCAN

 //**************************************************************
void ANALCANos(double usigm, double& usigm1) {
	//усиление каналов-400
	//double uv_f,u_f;
	//  int m=7;//razrjadnost
	int m_f;
	int nomkontip = 2;          //  Nc=2;
								//uv_f=300.*1.25*Ksl[0]*usigm*1.25; //150=355
								//UKON(55,nomkontip,KT,KN,usigm,usigm1);uv_f=usigm1;
								//	uv_f=.625*UUOS( uv_f)*1;
	uv_f = .625 * 355 * usigm;
	for (m_f = 0; m_f < 3; m_f++) {
		UKON(30 + m_f, nomkontip, KT, KN, uv_f, usigm1);
		uv_f = usigm1;
	}//m_f
	 //cout<<"n="<<n<<"  m="<<m<<"\n";
	 //getch();
}//ANALCANos

 //  --------------------------------
 /*
 // 8 mV
 u_acp1[0][1]=0.008041;  u_acp2[0][1]=0.007935;  u_acp1[1][1]=0.007904;  u_acp2[1][1]=0.007187;
 u_acp1[0][2]=0.007248;  u_acp2[0][2]=0.00766;  u_acp1[1][2]=0.006989;  u_acp2[1][2]=0.006821;
 u_acp1[0][3]=0.007446;  u_acp2[0][3]=0.007507;  u_acp1[1][3]=0.007126;  u_acp2[1][3]=0.006989;
 u_acp1[0][4]=0.008087;  u_acp2[0][4]=0.007584;  u_acp1[1][4]=0.008026;  u_acp2[1][4]=0.00737;
 u_acp1[0][5]=0.007324;  u_acp2[0][5]=0.007782;  u_acp1[1][5]=0.007813;  u_acp2[1][5]=0.00737;
 u_acp1[0][6]=0.007401;  u_acp2[0][6]=0.007584;  u_acp1[1][6]=0.007019;  u_acp2[1][6]=0.007492;
 u_acp1[0][7]=0.007172;  u_acp2[0][7]=0.007614;  u_acp1[1][7]=0.007523;  u_acp2[1][7]=0.00769;
 u_acp1[0][8]=0.007294;  u_acp2[0][8]=0.007477;  u_acp1[1][8]=0.00737;  u_acp2[1][8]=0.008026;
 u_acp1[0][9]=0.007736;  u_acp2[0][9]=0.006989;  u_acp1[1][9]=0.007278;  u_acp2[1][9]=0.007614;
 u_acp1[0][10]=0.007675;  u_acp2[0][10]=0.007126;  u_acp1[1][10]=0.007965;  u_acp2[1][10]=0.007187;
 //	      */
 //-----------------------------------------------------------
 /*
 //0 mV
 u_acp1[0][1]=0.00058;  u_acp2[0][1]=0.000259;  u_acp1[1][1]=0.000427;  u_acp2[1][1]=-0.00029;
 u_acp1[0][2]=-0.000259;  u_acp2[0][2]=0.000153;  u_acp1[1][2]=-0.00061;  u_acp2[1][2]=-0.000671;
 u_acp1[0][3]=-9.155273e-05;  u_acp2[0][3]=-0.000122;  u_acp1[1][3]=-0.000427;  u_acp2[1][3]=-0.000458;
 u_acp1[0][4]=0.00058;  u_acp2[0][4]=6.103516e-05;  u_acp1[1][4]=0.000488;  u_acp2[1][4]=-0.000168;
 u_acp1[0][5]=-0.000153;  u_acp2[0][5]=0.00032;  u_acp1[1][5]=0.00029;  u_acp2[1][5]=-0.000107;
 u_acp1[0][6]=-0.000259;  u_acp2[0][6]=0.000137;  u_acp1[1][6]=-0.000473;  u_acp2[1][6]=-4.577637e-05;
 u_acp1[0][7]=-0.000336;  u_acp2[0][7]=0.000107;  u_acp1[1][7]=7.629395e-05;  u_acp2[1][7]=0.000168;
 u_acp1[0][8]=-0.000183;  u_acp2[0][8]=-1.525879e-05;  u_acp1[1][8]=-9.155273e-05;  u_acp2[1][8]=0.000565;
 u_acp1[0][9]=0.000198;  u_acp2[0][9]=-0.000534;  u_acp1[1][9]=-0.000244;  u_acp2[1][9]=0.000168;
 u_acp1[0][10]=0.000153;  u_acp2[0][10]=-0.000397;  u_acp1[1][10]=0.000412;  u_acp2[1][10]=-0.000244;
 //		 */
 //--------------------------------
 /*
 //		 2 mV
 u_acp1[0][1]=0.002243;  u_acp2[0][1]=0.002075;  u_acp1[1][1]=0.002151;  u_acp2[1][1]=0.001648;
 u_acp1[0][2]=0.001694;  u_acp2[0][2]=0.001968;  u_acp1[1][2]=0.001648;  u_acp2[1][2]=0.001511;
 u_acp1[0][3]=0.001785;  u_acp2[0][3]=0.002014;  u_acp1[1][3]=0.001556;  u_acp2[1][3]=0.001434;
 u_acp1[0][4]=0.002365;  u_acp2[0][4]=0.002014;  u_acp1[1][4]=0.002197;  u_acp2[1][4]=0.001724;
 u_acp1[0][5]=0.001923;  u_acp2[0][5]=0.002075;  u_acp1[1][5]=0.001968;  u_acp2[1][5]=0.001755;
 u_acp1[0][6]=0.001678;  u_acp2[0][6]=0.001938;  u_acp1[1][6]=0.001587;  u_acp2[1][6]=0.001755;
 u_acp1[0][7]=0.001709;  u_acp2[0][7]=0.001877;  u_acp1[1][7]=0.001694;  u_acp2[1][7]=0.001907;
 u_acp1[0][8]=0.001587;  u_acp2[0][8]=0.001877;  u_acp1[1][8]=0.001831;  u_acp2[1][8]=0.002335;
 u_acp1[0][9]=0.001785;  u_acp2[0][9]=0.001343;  u_acp1[1][9]=0.001801;  u_acp2[1][9]=0.002136;
 u_acp1[0][10]=0.001663;  u_acp2[0][10]=0.001862;  u_acp1[1][10]=0.002274;  u_acp2[1][10]=0.001709;
 //           */
 //--------------------------------
 /*
 //		4 mV
 u_acp1[0][1]=0.00386;  u_acp2[0][1]=0.004074;  u_acp1[1][1]=0.004333;  u_acp2[1][1]=0.003464;
 u_acp1[0][2]=0.003601;  u_acp2[0][2]=0.003815;  u_acp1[1][2]=0.003586;  u_acp2[1][2]=0.003433;
 u_acp1[0][3]=0.003601;  u_acp2[0][3]=0.003723;  u_acp1[1][3]=0.003357;  u_acp2[1][3]=0.003357;
 u_acp1[0][4]=0.00441;  u_acp2[0][4]=0.003998;  u_acp1[1][4]=0.003922;  u_acp2[1][4]=0.00354;
 u_acp1[0][5]=0.003799;  u_acp2[0][5]=0.004105;  u_acp1[1][5]=0.003845;  u_acp2[1][5]=0.003708;
 u_acp1[0][6]=0.003235;  u_acp2[0][6]=0.004105;  u_acp1[1][6]=0.003464;  u_acp2[1][6]=0.00386;
 u_acp1[0][7]=0.003372;  u_acp2[0][7]=0.003967;  u_acp1[1][7]=0.003571;  u_acp2[1][7]=0.003769;
 u_acp1[0][8]=0.003494;  u_acp2[0][8]=0.003799;  u_acp1[1][8]=0.003662;  u_acp2[1][8]=0.004105;
 u_acp1[0][9]=0.003708;  u_acp2[0][9]=0.003357;  u_acp1[1][9]=0.003616;  u_acp2[1][9]=0.003738;
 u_acp1[0][10]=0.003418;  u_acp2[0][10]=0.003662;  u_acp1[1][10]=0.003983;  u_acp2[1][10]=0.003693;
 //         */
 //---------------------------------------------------
 /*
 //		   6 mV
 u_acp1[0][1]=0.005981;  u_acp2[0][1]=0.00589;  u_acp1[1][1]=0.006027;  u_acp2[1][1]=0.005402;
 u_acp1[0][2]=0.005219;  u_acp2[0][2]=0.005875;  u_acp1[1][2]=0.005203;  u_acp2[1][2]=0.005051;
 u_acp1[0][3]=0.005493;  u_acp2[0][3]=0.005569;  u_acp1[1][3]=0.005234;  u_acp2[1][3]=0.005264;
 u_acp1[0][4]=0.006134;  u_acp2[0][4]=0.005798;  u_acp1[1][4]=0.005981;  u_acp2[1][4]=0.005341;
 u_acp1[0][5]=0.005737;  u_acp2[0][5]=0.005875;  u_acp1[1][5]=0.005692;  u_acp2[1][5]=0.005585;
 u_acp1[0][6]=0.005219;  u_acp2[0][6]=0.005676;  u_acp1[1][6]=0.005463;  u_acp2[1][6]=0.005844;
 u_acp1[0][7]=0.005417;  u_acp2[0][7]=0.005859;  u_acp1[1][7]=0.005569;  u_acp2[1][7]=0.005722;
 u_acp1[0][8]=0.005447;  u_acp2[0][8]=0.005692;  u_acp1[1][8]=0.00563;  u_acp2[1][8]=0.006226;
 u_acp1[0][9]=0.00589;  u_acp2[0][9]=0.005096;  u_acp1[1][9]=0.005447;  u_acp2[1][9]=0.005753;
 u_acp1[0][10]=0.005661;  u_acp2[0][10]=0.005463;  u_acp1[1][10]=0.006165;  u_acp2[1][10]=0.005341;
 //           */
 //--------------------------------
 /*
 //		    12 mV
 u_acp1[0][1]=0.011414;  u_acp2[0][1]=0.011658;  u_acp1[1][1]=0.011887;  u_acp2[1][1]=0.011002;
 u_acp1[0][2]=0.011063;  u_acp2[0][2]=0.011398;  u_acp1[1][2]=0.011063;  u_acp2[1][2]=0.011002;
 u_acp1[0][3]=0.011139;  u_acp2[0][3]=0.01123;  u_acp1[1][3]=0.01088;  u_acp2[1][3]=0.010849;
 u_acp1[0][4]=0.011856;  u_acp2[0][4]=0.011505;  u_acp1[1][4]=0.011459;  u_acp2[1][4]=0.011032;
 u_acp1[0][5]=0.011307;  u_acp2[0][5]=0.011627;  u_acp1[1][5]=0.011276;  u_acp2[1][5]=0.01123;
 u_acp1[0][6]=0.010727;  u_acp2[0][6]=0.011597;  u_acp1[1][6]=0.010986;  u_acp2[1][6]=0.011414;
 u_acp1[0][7]=0.010956;  u_acp2[0][7]=0.011475;  u_acp1[1][7]=0.011108;  u_acp2[1][7]=0.011215;
 u_acp1[0][8]=0.011124;  u_acp2[0][8]=0.011307;  u_acp1[1][8]=0.011124;  u_acp2[1][8]=0.011658;
 u_acp1[0][9]=0.011307;  u_acp2[0][9]=0.01091;  u_acp1[1][9]=0.011093;  u_acp2[1][9]=0.011292;
 u_acp1[0][10]=0.010941;  u_acp2[0][10]=0.011154;  u_acp1[1][10]=0.011612;  u_acp2[1][10]=0.011169;
 //       */
 //--------------------------------
 /*
 //		     16 mV
 u_acp1[0][1]=0.015594;  u_acp2[0][1]=0.01535;  u_acp1[1][1]=0.015427;  u_acp2[1][1]=0.014755;
 u_acp1[0][2]=0.014755;  u_acp2[0][2]=0.015152;  u_acp1[1][2]=0.014511;  u_acp2[1][2]=0.014389;
 u_acp1[0][3]=0.014969;  u_acp2[0][3]=0.015106;  u_acp1[1][3]=0.01474;  u_acp2[1][3]=0.014526;
 u_acp1[0][4]=0.015579;  u_acp2[0][4]=0.015045;  u_acp1[1][4]=0.015503;  u_acp2[1][4]=0.014908;
 u_acp1[0][5]=0.014969;  u_acp2[0][5]=0.015228;  u_acp1[1][5]=0.015289;  u_acp2[1][5]=0.014938;
 u_acp1[0][6]=0.014862;  u_acp2[0][6]=0.015106;  u_acp1[1][6]=0.014603;  u_acp2[1][6]=0.015045;
 u_acp1[0][7]=0.014771;  u_acp2[0][7]=0.015091;  u_acp1[1][7]=0.015076;  u_acp2[1][7]=0.015152;
 u_acp1[0][8]=0.014771;  u_acp2[0][8]=0.015045;  u_acp1[1][8]=0.014954;  u_acp2[1][8]=0.015549;
 u_acp1[0][9]=0.015244;  u_acp2[0][9]=0.014511;  u_acp1[1][9]=0.014832;  u_acp2[1][9]=0.015091;
 u_acp1[0][10]=0.015244;  u_acp2[0][10]=0.014587;  u_acp1[1][10]=0.015533;  u_acp2[1][10]=0.014694;
 //            */
 //---------------------------------------------------
 /*
 //		       20 mV
 u_acp1[0][1]=0.018982;  u_acp2[0][1]=0.01915;  u_acp1[1][1]=0.019424;  u_acp2[1][1]=0.018524;
 u_acp1[0][2]=0.018524;  u_acp2[0][2]=0.018951;  u_acp1[1][2]=0.018555;  u_acp2[1][2]=0.018448;
 u_acp1[0][3]=0.018646;  u_acp2[0][3]=0.018814;  u_acp1[1][3]=0.018417;  u_acp2[1][3]=0.018265;
 u_acp1[0][4]=0.019333;  u_acp2[0][4]=0.018982;  u_acp1[1][4]=0.018967;  u_acp2[1][4]=0.018555;
 u_acp1[0][5]=0.018875;  u_acp2[0][5]=0.01918;  u_acp1[1][5]=0.018784;  u_acp2[1][5]=0.018768;
 u_acp1[0][6]=0.01825;  u_acp2[0][6]=0.019043;  u_acp1[1][6]=0.0186;  u_acp2[1][6]=0.018906;
 u_acp1[0][7]=0.018402;  u_acp2[0][7]=0.018982;  u_acp1[1][7]=0.018646;  u_acp2[1][7]=0.018753;
 u_acp1[0][8]=0.018677;  u_acp2[0][8]=0.018768;  u_acp1[1][8]=0.018661;  u_acp2[1][8]=0.019333;
 u_acp1[0][9]=0.018982;  u_acp2[0][9]=0.018372;  u_acp1[1][9]=0.0186;  u_acp2[1][9]=0.01889;
 u_acp1[0][10]=0.018478;  u_acp2[0][10]=0.0186;  u_acp1[1][10]=0.019196;  u_acp2[1][10]=0.018692;
 //			*/
 //--------------------------------
 /*
 //	   30 mV
 u_acp1[0][1]=0.028717;  u_acp2[0][1]=0.028427;  u_acp1[1][1]=0.028641;  u_acp2[1][1]=0.028;
 u_acp1[0][2]=0.027924;  u_acp2[0][2]=0.028519;  u_acp1[1][2]=0.027588;  u_acp2[1][2]=0.027588;
 u_acp1[0][3]=0.028046;  u_acp2[0][3]=0.028107;  u_acp1[1][3]=0.027725;  u_acp2[1][3]=0.027924;
 u_acp1[0][4]=0.028656;  u_acp2[0][4]=0.028397;  u_acp1[1][4]=0.028717;  u_acp2[1][4]=0.027954;
 u_acp1[0][5]=0.028198;  u_acp2[0][5]=0.02832;  u_acp1[1][5]=0.02832;  u_acp2[1][5]=0.028122;
 u_acp1[0][6]=0.027847;  u_acp2[0][6]=0.028214;  u_acp1[1][6]=0.028076;  u_acp2[1][6]=0.028336;
 u_acp1[0][7]=0.027786;  u_acp2[0][7]=0.028412;  u_acp1[1][7]=0.028214;  u_acp2[1][7]=0.028259;
 u_acp1[0][8]=0.02803;  u_acp2[0][8]=0.028259;  u_acp1[1][8]=0.028275;  u_acp2[1][8]=0.028687;
 u_acp1[0][9]=0.028381;  u_acp2[0][9]=0.027664;  u_acp1[1][9]=0.027954;  u_acp2[1][9]=0.028183;
 u_acp1[0][10]=0.028198;  u_acp2[0][10]=0.028015;  u_acp1[1][10]=0.028656;  u_acp2[1][10]=0.027939;
 */