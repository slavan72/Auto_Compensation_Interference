# if !defined(__F4D_LIB_H)
#define __F4D_LIB_H

//#include <math.h>
//********************************************************
void TEFI(float umo, float azo, float um, float az, float&tet, float&fi){ 
	//Переход из системы (УГОЛ МЕСТА,АЗИМУТХ[град]) в сферическую
   //связанную с полотном

   //umo и azo - ум и аз нормали к полотну(град)
   //uv и az  -  ум и аз точки наблюдения(град)
   //tet и fi -  сферич.координаты(град) точки наблюдения в сист.коорд.
   //            полотна(ось Z-нормаль к полотну,ось X -горизонтальна)
   float x, y, z, sa, ca, su, cu, suo, cuo, gr; //double
   gr=M_PI/180.;
   sa=(az-azo)*gr;  
   ca=cos(sa); 
   sa=sin(sa);
   su=um*gr;  
   cu=cos(su); 
   su=sin(su);
   suo=umo*gr;  
   cuo=cos(suo); 
   suo=sin(suo);
   x=cu*sa;
   y=su*cuo-cu*suo*ca;
   z=su*suo+cu*cuo*ca;
   tet=acos(z)/gr;
   fi=ARCTG(x,y)/gr;
  }//TEFI
//********************************************************
//             БОЛЬШОЕ ПОЛОТНО
//-----------------------------------------------------------
	struct geom_polotno{
		float dx;
		float dy;
		float gek;
	}geom_40;
	//-----------------------------------------------------------
	int const TIPFAR=3;   //количество типов реш
	//float G[TIPFAR]={0.,-15.5,-11.};//усиление  типа
	//float G[TIPFAR]={0.,-17.,-13.};//усиление  типа (основной вариант)
	//float G[TIPFAR]={0.,-17.5,-17.5};//усиление  типа
	//float G[TIPFAR]={0.,-20.5,-14.};//усиление  типа
	float G[TIPFAR]={0.,-17.,-14.};//усиление  типа (Vorobjev)
	float normir[TIPFAR];
	int ampn[TIPFAR];
	unsigned char (*amp[TIPFAR]);
	//-----------------------------------------------------------
	unsigned char I_40=112; //строки
	unsigned char J_40[112]=
				  {64,64,64,64,64,64,64,64, 80,80,80,80,80,80,80,80,
				   96,96,96,96,96,96,96,96, 96,96,96,96,96,96,96,96,
				   96,96,96,96,96,96,96,96, 96,96,96,96,96,96,96,96,
				   96,96,96,96,96,96,96,96, 96,96,96,96,96,96,96,96,
				   96,96,96,96,96,96,96,96, 96,96,96,96,96,96,96,96,
				   96,96,96,96,96,96,96,96, 96,96,96,96,96,96,96,96,
				   80,80,80,80,80,80,80,80, 64,64,64,64,64,64,64,64 };


	unsigned int Nelem_40;

void NELEM_F4(unsigned char I_40,unsigned char J_40[], unsigned int& Nelem_40){ 
	unsigned char i;
   Nelem_40=0;
   for(i=0; i<I_40; i++)
     {
		 Nelem_40=Nelem_40+J_40[i];
       }//i
  }// NELEM_F4
//-----------------------------------------------------------

// void AMP0(unsigned char *amp0,float&normir)
 void AMP0(unsigned char *amp[], float& normir){
	 float e, h, st, ct, sf, cf, st0, sf0, cf0;         
	 int i, j, m, M;
     float x, y, z, r, ro, a, b, gek, u, v, lam;                  
     float aaa;
//112 строк,96 столбцов,модуль горизонтально
       lam=300./9.5;

		 //double const pi=3.1415926;  
		 double const k=2.*M_PI;
//         float dx=.64917;float dy=.53833; gek=.25; корея
		 float dx=.75; 
		 float dy=.666; 
		 gek=.25; //f400

	//  80        55         1800 корея
	//  90        58         3680 ММ
		a=90./lam;b=70./lam; 
		ro=3150./lam;      //     ro=3100./lam;
		y=-dy*((float)I_40/2.+.5);    
		e=M_PI*M_PI/8.*ro*ro;
		normir=0.;
		m=-1;

		for(i=0; i<I_40; i++){
			x=-dx*((float)J_40[i]/2.+.5)-dx*gek;gek=-gek; 
			y=y+dy;

			  for(j=0; j<J_40[i]; j++){
				  m=m+1;
				//m=i+J*j;
				x=x+dx;
				z=x*x+y*y;
				r=sqrt(ro*ro+z);
				u=k*y/r*a/2.; 
				v=k*x/r*b/2.;
				u=cos(u)/(u*u-.25*M_PI*M_PI); 
				v=sin(v)/v; 
				ct=ro/r; // ct=1./ct/ct/sqrt(sqrt(ct));
					  // aaa=(x*x+y*y*ct)/z;
				h=-(1.+ct)*u*v/r*(x*x+y*y*ct)/(z*r) *e;
				  *(amp[0]+m)=(unsigned char)floor(255.*h+.5);
				normir=normir+ *(amp[0]+m);
			 }  //i
	    }  //j
    } //AMP0
//-----------------------------------------------------------

void DA0(float t0, float f0, float t, float f, float normir,
//	 double&re,double&im)
	 unsigned char *amp[], double&re, double&im){

	float c, c1, e, e1, h, st, ct, ct3, sf, cf, st0, sf0, cf0;         
	int i, j, m;
    //double const pi=3.1415926; 
	double const k=2.*M_PI;
//     float dx=.53833;float dy=.64917; float gek=.25;//корея
    float dx=.75;
	float dy=.6666; 
	float gek=.25;//f400
    float x, y, z, r, r1, ro, r2, a, b, u, v, sig;   
	double a1;
   //double an,an1;
   //an=3579480.1*(3+t0);   an1=46823373.1*(6.5+f0);
    sig=.3;
	sig=.0;
//double are=1e8*fabs(4.37*t0+6.62*f0);
    double aim=1e8*fabs(5.12*t0+3.91*f0+3.45287852);
// double fon=.00001/1.8;//уровень фона
// int N=6784;  //количество элементов
// sig=sqrt(3.*log(1.+N*fon));sig=0.53; float dg=exp(sig*sig/6.)*1.004;
    float dg;//dg=4.*sig*sig/12.; dg=exp(dg/2.)*1.01041;
// dg=1./(1.-exp(-sig*sig));
	dg=exp(sig*sig*.5);
	st=sin(t);ct=cos(t);  
	ct3=sqrt(ct*ct*ct); // ct3=1;
	sf=sin(f);
	cf=cos(f);
	st0=sin(t0);
	sf0=sin(f0);
	cf0=cos(f0);
	y=-dy*(float)(I_40/2+.5);    
	e=c=e1=c1=0.;
	m=-1;
	r1=0;
	    //NRV1(an);
	    RND(aim);
	for(i=0; i<I_40; i++){
		x=-dx*(float)(J_40[i]/2+.5)-dx*gek;
		gek=-gek; 
		y=y+dy;
	  for(j=0;j<J_40[i];j++){
		  x=x+dx;
		  r=-k*(x*(st*cf-st0*cf0)+y*(st*sf-st0*sf0));//+sig*NRV1(an);//RND(an);
		  r2=-k*(x*st*cf+y*st*sf+1000*RND(aim));//+sig*NRV1(an);//RND(an);
			//r2=1000*k*RND(aim);//+1000*r;
		  m=m+1;
		//m=i+j*J;
		//h=amp0[m];
		  e=e+*(amp[0]+m)*sin(r);
		  c=c+*(amp[0]+m)*cos(r);
		  e1=e1+sin(r2);  c1=c1+cos(r2);
	//cout<<m<<"  "<<r2<<"\n";getch();
		r1=r1+r;
	  }  //i
    }  //j
//re=(float)(c*dg/normir*ct3+.0005*NRV1(1e8*fabs(cos(r1))));
//im=(float)(e*dg/normir*ct3+.0005*NRV1(1e8*fabs(sin(r1))));
	re=(float)(c*dg/normir*ct3*1+.000051*c1);
	im=(float)(e*dg/normir*ct3*1+.000051*e1);
//cout<<"           "<<re*re+im*im<<"\n";getch();
}  //конец DA0

//**************************************************************************
//*************************************************************************
//                ПОЛОТНО 16*16
//--------------------------------------------------------------------------
//-----------------------------------------------------------
unsigned char I_41=16; //строки
unsigned char J_41[16]={16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16};
unsigned int Nelem_41;
//-----------------------------------------------------------
//              амплитудное распределение
//-----------------------------------------------------------

 void AMP1(unsigned char *amp[], float&normir){
	 float h, st, ct, sf, cf, st0, sf0, cf0;         
	 int i, j, m, M;
	 float x, y, z, r, ro, a, b, gek, u, v, lam;                  
	 float aaa;
//112 строк,96 столбцов,модуль горизонтально
     lam=300./9.5;
// double const pi=3.1415926; //double const k=2.*pi;
//     float dx=.64917;float dy=.53833; gek=.25; корея
     float dx=22.5/lam;
	 float dy=20./lam; 
	  gek=.25; //f400
 //  80        55         1800 корея
 //  90        58         3680 ММ
//a=90./lam;b=70./lam; ro=3150./lam;      //     ro=3100./lam;
  // float  rx=5.9;
  // float  ry=5.;
	ro=7.5;
	y=-dy*((float)I_41/2.+.5);   // e=pi*pi/8.*ro*ro;
	normir=0.;
		m=-1;

	for(i=0; i<I_41; i++){
		x=-dx*((float)J_41[i]/2.+.5)-dx*gek;
		gek=-gek;
		y=y+dy;
	  for(j=0;j<J_41[i];j++){
			m=m+1;
			//m=i+J*j;
			x=x+dx;
			z=sqrt(x*x+y*y*1.27);
			//h=(.15+(1.+cos(1.*pi*z/ro)))/(.15+2.);
			h=.28+pow((1.-pow(fabs(z/ro),2)),2);//+.29*(j-8)/8.;
			   h=h/1.28;
			*(amp[1]+m)=(unsigned char)(255*h+.5);
			normir=normir+ *(amp[1]+m);
		}  //i
	}  //j
} //AMP1

//-----------------------------------------------------------
void DA1(float t0, float f0, float t, float f, float normir,
	 unsigned char *amp[], double&re, double&im){
		 float c, e, h, st, ct, ct3, sf, cf, st0, sf0, cf0, lam;         
		 int i, j, m;
 //double const pi=3.1415926; 
 //double const k=2.*M_PI;
	       lam=300./9.5;
//     float dx=.53833;float dy=.64917; float gek=.25;//корея
         float dx=22.5/lam;float dy=20./lam; 
		 float gek=.25;//f400
		 float x, y, z, r, ro, a, b, u, v, sig;  
		 double an;
		 double a1;
		 an=35794801;
		 sig=.5;                                   //    sig=0;
// double fon=.00001/1.8;//уровень фона
// int N=6784;  //количество элементов
// sig=sqrt(3.*log(1.+N*fon));sig=0.53; float dg=exp(sig*sig/6.)*1.004;
		 float dg;//dg=4.*sig*sig/12.; dg=exp(dg/2.)*1.01041;
// dg=1./(1.-exp(-sig*sig));
		 dg=exp(sig*sig*.5);
		 st=sin(t);
		 ct=cos(t);    
		 ct3=sqrt(ct*ct*ct);
		 sf=sin(f);
		 cf=cos(f);
		 st0=sin(t0);
		 sf0=sin(f0);
		 cf0=cos(f0);
		 y=-dy*(float)(I_41/2+.5);    
		 e=c=0.;
		 m=-1;

	for(i=0; i<I_41; i++){
		x=-dx*(float)(J_41[i]/2+.5)-dx*gek; 
		gek=-gek; 
		y=y+dy;
	  for(j=0;j<J_41[i];j++){
		   x=x+dx; 
		   r=-2.*M_PI*(x*(st*cf-st0*cf0)+y*(st*sf-st0*sf0))+sig*NRV1(an);//RND(an);
		   m=m+1;
		//m=i+j*J;
		//h=amp0[m];
		   e=e+*(amp[1]+m)*sin(r);
		   c=c+*(amp[1]+m)*cos(r);
		}  //i
	}  //j
	re=(float)(c*dg/normir*ct3); 
	im=(float)(e*dg/normir*ct3);
}  //конец DA0

//*************************************************************************
//*************************************************************************
//                ПОЛОТНО 48*16
//--------------------------------------------------------------------------
//-----------------------------------------------------------
unsigned char I_42=16; //строки
unsigned int Jx42=48; //столбцы
unsigned char J_42[16]=//{16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16};
		       {48,48,48,48,48,48,48,48,48,48,48,48,48,48,48,48};
// for(j=0;j<Jx42;j++)J_42[j]=Jx42;
unsigned int Nelem_42;
//-----------------------------------------------------------
//              амплитудное распределение
//-----------------------------------------------------------

 void AMP2(unsigned char *amp[],float&normir){
	 float e,h,st,ct,sf,cf,st0,sf0,cf0;         
	 int i, j, M; // m,
	 int m = 0;
	 float x, y, z, r, ro, a, b, gek, u, v, lam;                  
	 float aaa;
//112 строк,96 столбцов,модуль горизонтально
     lam=300./9.5;
// double const pi=3.1415926; //double const k=2.*pi;
//     float dx=.64917;float dy=.53833; gek=.25; корея
      float dx=22.5/lam;float dy=20./lam; gek=.25; //f400
 //  80        55         1800 корея
 //  90        58         3680 ММ
//a=90./lam;b=70./lam; ro=3150./lam;      //     ro=3100./lam;
     float  ax=1;
     float  ay=1.;
     ro=0;
     y=-dy*((float)I_42/2.+.5);

	for(i=0; i<I_42; i++){
		x=-dx*((float)J_42[i]/2.+.5)-dx*gek;
		gek=-gek; 
		y=y+dy;
	  for(j=0; j<J_42[i];j ++){
		  m=m+1;
		//m=i+J*j;
		 x=x+dx;
		 z=sqrt(x*x*1+y*y*ay);
		  if(z>ro) ro=z;
	  } //j
   } //i

	y=-dy*((float)I_42/2.+.5);
	normir=0.;
	m=-1;

	for(i=0; i<I_42; i++) {
		x=-dx*((float)J_42[i]/2.+.5)-dx*gek;
		gek=-gek; 
		y=y+dy;
	  for(j=0;j<J_42[i];j++){
		  m=m+1;
		//m=i+J*j;
		  x=x+dx;
		  z=sqrt(x*x*ax+y*y*ay);
		//h=(.15+(1.+cos(1.*pi*z/ro)))/(.15+2.);
	   // h=.15+pow(fabs(1.-pow(fabs(z/ro),1.8)),2);//+.29*(j-8)/8.;
	   //    h=h/1.15;
		  h=.1+exp(-x*x/289*4.-y*y/40*3.);      
		  h=h/1.1;
		  *(amp[2]+m)=(unsigned char)(255*h+.5);
		  normir=normir+ *(amp[2]+m);
	  }  //j
    }  //i
 } //AMP2

//-----------------------------------------------------------
void DA2(float t0,float f0,float t,float f,float normir,
	 unsigned char *amp[],double&re,double&im){
		 float c, e, h, st, ct, ct3, sf, cf, st0, sf0, cf0, lam;         
		 int i, j, m;
		 //double const pi=3.1415926; 
		 //double const k=2.*M_PI;
	       lam=300./9.5;
//     float dx=.53833;float dy=.64917; float gek=.25;//корея
		float dx=22.5/lam;
		float dy=20./lam; 
		float gek=.25;//f400
		float x,y,z,r,ro,a,b,u,v,sig;  
		double an; 
		double a1;
		an=35794801;
		sig=.5;               //                        sig=0;
// double fon=.00001/1.8;//уровень фона
		float dg;//dg=4.*sig*sig/12.; dg=exp(dg/2.)*1.01041;
// dg=1./(1.-exp(-sig*sig));
		dg=exp(sig*sig*.5);
		st=sin(t); ct=cos(t);   ct3=sqrt(ct*ct*ct);
		sf=sin(f);cf=cos(f);st0=sin(t0);sf0=sin(f0);cf0=cos(f0);
		y=-dy*(float)(I_42/2+.5);    e=c=0.;
		m=-1;

	for(i=0; i<I_42; i++){
		x=-dx*(float)(J_42[i]/2+.5)-dx*gek;gek=-gek; 
		y=y+dy;
	  for(j=0; j<J_42[i]; j++){
		  x=x+dx; 
		  r=-2.0*M_PI*(x*(st*cf-st0*cf0)+y*(st*sf-st0*sf0))+sig*NRV1(an);//RND(an);
		  m=m+1;
		//m=i+j*J;
		//h=amp0[m];
		  e=e+*(amp[2]+m)*sin(r);
		  c=c+*(amp[2]+m)*cos(r);
	   }  //i
	}  //j
	re=(float)(c*dg/normir*ct3); im=(float)(e*dg/normir*ct3);
}  //конец DA2

//*************************************************************************
//*************************************************************************
//      расчет амп распр поля для всеж типов реш

 void AMP400(int tipfar,unsigned char *amp[],float normir[]){ 
	if (tipfar==0) AMP0( amp,normir[tipfar]); //tipfar=0

    if (tipfar==1) AMP1( amp,normir[tipfar]); //tipfar=1

    if (tipfar==2) AMP2( amp,normir[tipfar]); //tipfar=2

    }//AMP
//**************************************************************************
//      расчет ДН для всеж типов реш
void D400(float t0,float f0,float t,float f,float normir[],
	  int tipfar,unsigned char *amp[],double&re,double&im){
  if(tipfar==0) DA0(t0,f0,t,f,normir[tipfar],amp,re,im);

  if(tipfar==1) DA1(t0,f0,t,f,normir[tipfar],amp,re,im);

  if(tipfar==2) DA2(t0,f0,t,f,normir[tipfar],amp,re,im);
 }//D400

//********************************************************
			 /*
 void ZADERG(int nom,int nzad,long int n,double u1,double&u2)
 {if(nzad>0){u2=zader[nom][n%(nzad)];zader[nom][n%(nzad)]=u1;}
   else u2=u1;
  }// ZADERG
			   */
 float fmin0=52.5e6; float fmax0=67.5e6; float fo0;//полоса помех на входе
 float fmin1=57.5e6; float fmax1=62.5e6; float fo1;//полоса рез ус
 float fmin2=57.5e6; float fmax2=62.5e6; float fo2;//полоса рез ус OC
// float fmin3=51.5e6;float fmax3=70e6;float fo3;//полоса REGUS
// float fmin3=53e6;float fmax3=68e6;float fo3;//полоса REGUS    29.8 db
// float fmin3=54.5e6;float fmax3=66.0e6;float fo3;//полоса REGUS   27 db
 float fmin3=53.35e6; float fmax3=67.3e6; float fo3;//полоса REGUS   27 db

 float fmin4=55.5e6; float fmax4=59.8e6; float fo4;//полоса OKFIL
 float fmin5=60.2e6; float fmax5=64.6e6; float fo5;//полоса OKFIL
 float fmin6=60.34e6; float fmax6=62.7e6; float fo6;//полоса OKFIL

 float fmin7=51e6; float fmax7=59e6; float fo7;//полоса VHODFIL_2
 float fmin8=61e6; float fmax8=70e6; float fo8;//полоса VHODFIL_2
 float fmin9=62.e6; float fmax9=66e6; float fo9;//полоса VHODFIL_2

//--------------------------------------------------------------------------
void REALIZ(int nr,int Np,int Nc,double v[4][5][2],double af[4 ][5][2],
	    double vs[5][2],double afs[5][2],double tns[5] ){
	double b,c;
    c=1e9+1e2*nr;
    int j,i,m,q;

  for(j=0; j<Np; j++)
	  for(i=0; i<Nc; i++)
		  af[j][i][0]=af[j][i][1]=0.0;

  for(j=0;j<Np;j++){
    //	 c=1e9+1e2*nr;
    for(m=0;m<2;m++){
      b=1e9*(1+RND(c))+1e8;
      for(i=0;i<Nc;i++){
		  v[j][i][m]=b;
//cout<<"помеха="<<j<<"  канал="<<i<<"  комп="<<m<<"  v="<<v[j][i][m]<<"\n";getch();
	  }//i
    }//m
  }//j
//-------------------------------------------------
  c=1e9+50+100*nr;
  RND(c);

  for(i=0; i<Nc; i++){
	  b=RND(c);
      tns[i]=(1+b)*.5;
//cout<<"канал "<<i<<"  b="<<b<<"  'tн'="<<tns[i]<<"\n";   getch();
    }//i

  for(i=0;i<Nc;i++)
	  afs[i][0]=afs[i][1]=0.;

  for(i=0; i<Nc; i++){
	   for(m=0;m<2;m++){
				      // c=1e9+50+100*nr;
		   vs[i][m]=1e9*(1+RND(c))+1.03e8;
       }
//cout<<"канал "<<i<<"  "<<vs[i][0]<<"  "<<vs[i][1]<<"\n";   getch();
    }//i
  }//REALIZ

//--------------------------------------------------------------------------
//********************************************************
//****************************************************
	  //	       /*
//  ИНТЕГРИРУЮЩАЯ ЦЕПОЧКА
int const Ntipintc=2; //B REGUS: тип 0 и FD: тип 1
int const Nintc=32; //24-REGUS(тип 0); 8-FD(тип 1)

struct INTCTIP{
	double e[Ntipintc];
	float tau[Ntipintc];
}INTCT;

void INTCHT(float dt,float tau,int tip,INTCTIP&INTCT);
struct INTCNOM{
	double u0[Nintc];
	double u1[Nintc];
	double du[Nintc];
	double un[Nintc];
}INTCN;
//void UINTCn(INTCNOM&INTCN);
//void UINTC(int nomintc,int nomintctip,INTCTIP INTCT,INTCNOM&INTCN,
  //	    double uv,double&u);
	    //		       */
//****************************************************
//********************************************************************
void INTCHT(float dt,float tau,int tip,INTCTIP&INTCT){ 
	INTCT.tau[tip]=tau; //tau=3e-7[сек] -постоянная времени интегр.цепочки
    INTCT.e[tip]=exp(-dt/tau);
  }//INTC

void UINTCn(INTCNOM&INTCN){  
	int n;
     for(n=0; n<Nintc; n++)
      INTCN.du[n]=INTCN.u0[n]=INTCN.u1[n]=INTCN.un[n]=0.;
   } //INTCn

void UINTC(int nom,int tip,INTCTIP INTCT,INTCNOM&INTCN,double uv,double&u){
	INTCN.du[nom]=uv-INTCN.un[nom];
    INTCN.u0[nom]=INTCN.u0[nom]+INTCN.du[nom];
    INTCN.u1[nom]=(INTCN.u1[nom]-INTCN.du[nom])*INTCT.e[tip];
    u=INTCN.u0[nom]+INTCN.u1[nom];
    INTCN.un[nom]=uv;
  }//FLF

//**************************************************************
//****************************************************
//  ДИФФЕРЕНЦИРУЮЩАЯ ЦЕПОЧКА
int const Ntipdifc=1; //B REGUS: тип 1
int const Ndifc=8;
struct DIFCTIP{
	double e[Ntipdifc];
	float tau[Ntipdifc];
}DIFCT;
void DIFCHT(float dt, float tau, int tip, DIFCTIP&DIFCT);
struct DIFCNOM{
	double u0[Ndifc];
	double u1[Ndifc];
	double du[Ndifc];
	double un[Ndifc];
}DIFCN;

void UDIFCn(DIFCNOM&DIFCN);

void UDIFC(int nomdifc, int nomdifctip, DIFCTIP DIFCT, DIFCNOM&DIFCN, double uv, double&u);
//****************************************************
//**************************************************************
void SUMMATOR1(int Nc, double usignal[], double uws[], double uwc[], double&usigm){
	 int n;
	 usigm=usignal[0];
//  ZADERG(0,L,usignal[0],usigm);
//---------------------------------------------------
  for( n=1; n<Nc; n++){
     usigm=usigm+(uwc[n-1]+uws[n-1]);
     }//n
//usigm=usignal[0];
 }//SUMMATOR1
//*******************************************************
int const Ntipkon=10;
int const Nkon=56;

struct KONNOM{
	float y[Nkon];
	float y0[Nkon];
	float ys[Nkon];
	float yc[Nkon];
	double du[Nkon];
	double un[Nkon];
	      }KN;

struct KONTIP{
	float f0[Ntipkon];
	float Q[Ntipkon];
	float W[Ntipkon];
	float Ku[Ntipkon];
	float g[Ntipkon];
	double h0[Ntipkon];
	double hc[Ntipkon];
	double hs[Ntipkon];
	double ct[Ntipkon];
	double st[Ntipkon];
}KT;

void KHT(float dt,float f0,float Q,float W,float Ku,int nomkontip,KONTIP&KT);

void UKON(int nomkon, int nomkontip, KONTIP KT, KONNOM&KN, double uv, double& u);

//**************************************************************
void UKON(int nomkon,int nomkontip,KONTIP KT,KONNOM&KN,double uv,double& u){
  double af_1;
  KN.du[nomkon]=uv-KN.un[nomkon];
  KN.y0[nomkon]=KN.y0[nomkon]+KT.h0[nomkontip]*KN.du[nomkon];
  KN.y[nomkon]=KN.y0[nomkon]+KN.du[nomkon]*KT.hc[nomkontip]+
	       KN.yc[nomkon]*KT.ct[nomkontip]-KN.ys[nomkon]*KT.st[nomkontip];
  af_1=KN.yc[nomkon];
//  py=KN.yc[nomkon];
  KN.yc[nomkon]=KN.du[nomkon]*KT.hc[nomkontip]+KN.yc[nomkon]*KT.ct[nomkontip]-
	       KN.ys[nomkon]*KT.st[nomkontip];

  KN.ys[nomkon]=KN.du[nomkon]*KT.hs[nomkontip]+
		KN.ys[nomkon]*KT.ct[nomkontip]+af_1*KT.st[nomkontip];

  u=KN.y[nomkon]*KT.g[nomkontip]; //C[n]=y;    //напряжение на выходе контура
  KN.un[nomkon]=uv;
  }  //UKON

//******************************************************
//			     /*
//входной фильтр
void VHODFIL_2(int nom, int nc, KONTIP KT, KONNOM KN, double u1, double& u2){
	//nom - номер VHODFIL_2
   //однозвенный ФСИ, Ku=1
  int m_f, nomkon, s; //, uv_f
  double uv_f;
  int q=2; //количество контуров в фильтре
  int nomkontip=7; //(7,8,9)
  //nomkon=0; //(0-14)
  u2=0;     s=-1;
    for(m_f=0; m_f<q; m_f++){
		UKON(0+q*nc+m_f, nomkontip+m_f, KT, KN, u1, 
			uv_f);//???????????????????????
      s=-s;
      u2=u2+s*uv_f;
//cout<<"  usigm2="<<usigm2<<"  ";
      }//m
    u2=u2/.93;
  }//VHODFIL_2
    //			      */
//**************************************************************
//*************************************************************************
void UKONn(KONNOM&KN){ 
	int n;
    for(n=0; n<Nkon; n++)
    KN.y[n]=KN.y0[n]=KN.ys[n]=KN.yc[n]=KN.du[n]=KN.un[n]=0;
    }//UKONn

//*************************************************************************
//**************************************************************
//оконечный фильтр
void OKFIL(int nom,double u1,double& u2){ 
	//nom - номер OKFIL
   //3 связанных контура последовательно
  int q=6; //количество контуров в фильтре
  int nomkontip=4;
  int nomkon,s,m_f;
  double uv_f,u_f;
  nomkon=42; //(0-14)

  u2=0;     
  s=-1;
    for(m_f=0; m_f<2; m_f++){
	  UKON(nomkon+m_f+q*nom, nomkontip+m_f, KT, KN, u1, uv_f);
      s=-s;
      u2=u2+s*uv_f;
      }//m

  u_f=u2;
  u2=0;     
  s=-1;
    for(m_f=0; m_f<2; m_f++){
	  UKON(nomkon+2+m_f+q*nom, nomkontip+m_f, KT, KN, u_f, uv_f);
      s=-s;
      u2=u2+s*uv_f;
      }//m

  u_f=u2;
  u2=0;     
  s=-1;
    for(m_f=0;m_f<2;m_f++){
	  UKON(nomkon+4+m_f+q*nom,nomkontip+m_f,KT,KN,u_f,uv_f);
      s=-s;
      u2=u2+s*uv_f;
      }//m

  u2=u2/.87;
  }//OKFIL

//**************************************************************
//**************************************************************
//**************************************************************
void KALIBROVKA(float dt, float fmin, float fmax, int N, double&kalibr, int&Nsl, float&fo){
	//N-кол-во контуров формирующих входную полосу
  int m, i, j, q; 
  unsigned long int n;
  int nomkontip=0;
  double uv, u, u1, p, b, c;
  float Q;
  int M=7000; //кол-во изменений слу вел  1500
  int Nust;
  double v[2], af[2];
  v[0]=12648465; 
  v[1]=85674534;

  RALLEY(v,af);

  Nust=3*(int)floor(5./(fmax-fmin)/dt); //Nust*dt-время установки
  Nsl=(int)floor(.5/((fmax-fmin)*dt));//Nsl*dt- период обновления случ вел
//cout<<"калибровка:Nsl="<<Nsl<<"\n";
  M=M+Nust/Nsl+1;
  float W=1;  
  float Ku=1;
  UKONn(KN);
  p=b=c=0;
  n=0;

  for(m=0; m<M; m++){
	  RALLEY(v,af);
//if(m%10==0)cout<<"m="<<m<<"\n";
     for(i=0; i<Nsl; i++){
	   n=n+1;
       u1=uv=af[0]*sin(2*M_PI*fo*i*dt+af[1]);
//cout<<n<<" u1="<<u1<<"  ";getch();
     b=b+af[0]*af[0];
//-------------------------------------------------------------
  VHODFIL_2(0, 0, KT, KN, uv, u);
//cout<<n<<"  Uвх="<<uv<<"  Uвых="<<u<<"\n";getch();
    uv=u;

       OKFIL(0,u,uv);
//cout<<n<<"  "<<u<<"  uv="<<uv<<"\n";getch();
       if(n > Nust){
		   c=c+u1*u1;
		   p=p+uv*uv;
	   }
//cout<<" c="<<c<<" p="<<p<<"\n";getch();
       }//i
     }//m
    n=n-Nust;
//getch();
     kalibr=(float)sqrt(c/p);
     p=sqrt(p/n);
     b=sqrt(b/n);
     c=sqrt(c/n);
//cout<<"kalibr="<<kalibr<<"\nb="<<b<<"  c="<<c<<"  p="<<p<<"\n";
			  /*
cout<<"KALIBROVKA: "<<M<<" замеров\n"
      "тип контуров формируюших полосу входного сигнала- "<<nomkontip<<
      "\n    количество- "<<N<<" на канал\n"<<
      "полоса входного сигнала-  "<<fmin*1e-6<<"-"<<fmax*1e-6<<" МГц\n"<<
      "kalibr=" <<kalibr<<"  Nsl="<<Nsl<<
      "  Pср вх="<<c*c<<"  Pср вых="<<p*p<<"\n";
			  */
cout<<"KALIBROVKA: "<<M<<" замеров  " <<"  Nsl="<<Nsl<<"\n"
   << "полоса сигнала на входе блока- "<<fmin*1e-6<<"-"<<fmax*1e-6<<" МГц\n"
   <<"Pспектра="<<b*b<< "  Pср вх="<<c*c<<"  Pср вых="<<p*p<<"\n"
   <<"kalibr=" <<kalibr<<"\n";
  _getch();
  }//KALIBROVKA
//**************************************************************
#endif __F4D_LIB_H

