# if !defined(__WLIB4_H)
#define __WLIB4_H

//#include <math.h>
//**********************************************************
struct DUAMEL_drobi{
	int N;
	int J;
	double Ko;
	double Z[8]/*корни*/;
	double a[8];
    int sop[8];
} ;
	//FAZOSa;

//**********************************************************
double FACT(int n){
	//ФАКТОРИАЛ
  double a;
  int i;

  if(n==0||n==1){
	  a=1;
	  goto ZZ;
  }
  a=1;
  for(i=1;i<=n;i++) a=a*i;
  ZZ:
  return a;
  } //FACT

//**********************************************************
void PROIZ(int M, int N, double A[], int&J, double B[]){
	//ДИФФЕРЕНЦИРОВАНИЕ ПОЛИНОМА
  // M-порядок произв.
  // N-степень полинома
  // A-коэф-ты полинома
  // J-степень поли.производной
  // B-коэф-ты поли.производной
  int j;   
  J=N-M;

  for(j=0; j<=N; j++)  B[j]=0;
  for(j=0; j<=J; j++)  B[j]=A[j+M]*FACT(M+j)/FACT(j);

  }//PROIZ

//**********************************************************
void POLY(complex<double> x, int N, double A[], complex<double>&f){
	int n;
//cout<<"POLY:\n";getch();
  f=complex<double>(0,0);
  for(n=N;n>=0;n--){
	  f=f*x+A[n];
	//cout<<"n="<<n<<"  f="<<f<<"\n";
      }
  }//POLY

//**********************************************************
void DELEN2(int N, double A[], double u, double v, int&M, double B[]){ 
	int i,j; 
	double b[10];
//cout<<"вошли в DELEN2\n";getch();
   M=N-2;
   for(i=0;i<=N;i++)  b[i]=A[i];
//cout<<"mbn\n";
   for(i=N-2;i>=0;i--){
	   B[i]=b[i+2];
       b[i+1]=b[i+1]-u*B[i];
       b[i]=b[i]-v*B[i];
//cout<<"B["<<i<<"]="<<B[i]<<"\n";getch();
     }//i
  }//DELEN2

//**********************************************************
void DELEN1(int N, double A[], double v, int&M, double B[]){ 
	int n;
   M=N-1;
   B[N-1]=A[N];
   for(n=N-2; n>=0; n--) B[n]=A[n+1]+v*B[n+1];
   }//DELEN1

//**********************************************************
void BARSTOE(int N, double A[], complex<double>&x1, complex<double>&x2,int&so){
	double b[12],c[12],u,v,d,e,g;  
	int n,i;
	long int j;
  complex<double> f1,f2;
//cout<<"BARST\n";
//for(n=0;n<=N;n++)cout<<"n="<<n<<"  A="<<A[n]<<"\n"; getch();
  b[0]=b[1]=c[0]=c[1]=0.;
  u=v=0.;
  j=0;
  ZZ1:j=j+1;
//cout<<"hfsd 1\n";getch();
    for(i=2;i<=N+2;i++){ 
		b[i]=A[N-i+2]+u*b[i-1]+v*b[i-2];
       c[i]=b[i]+u*c[i-1]+v*c[i-2];
       }//i

//cout<<"hfsd 2\n";getch();
    d=(c[N]*c[N]-c[N+1]*c[N-1]);    
	d=1./d;
    u=u+(b[N+2]*c[N-1]-b[N+1]*c[N])*d;
    v=v+(b[N+1]*c[N+1]-b[N+2]*c[N])*d;
//cout<<"BARSTOE: j="<<j<<"  u="<<u<<"  v="<<v<<"\n";getch();
//cout<<"hfsd 3\n";getch();
    d=u*u/4.+v; 
	e=sqrt(fabs(d));
//cout<<"d="<<d<<"  e="<<e<<"\n";getch();
    if(d<0){
		x1=complex<double>(u/2,e);
		x2=complex<double>(u/2,-e);so=1;
	} else {
		x1=complex<double>(u/2-e,0.);
		x2=complex<double>(u/2+e,0.);
		so=0;
	}
//cout<<"hfsd 4\n";getch();
    POLY(x1,N,A,f1);    
	POLY(x2,N,A,f2);
//cout<<"hfsd 5\n";getch();
    e=sqrt(norm(f1)); 
	g=sqrt(norm(f2));
    d=fabs(1e-15*A[0]);
    if(j>100){
		cout<<"ПОНИЗЬ ТОЧНОСТЬ В 'BARSTOE'\n"; goto ZZ2;
	}
//cout<<"j="<<j<<"  x1="<<x1<<"  u="<<e/d<<"  x2="<<x2<<"  v="<<g/d<<"\n";getch();
    if(e>d||g>d)  goto ZZ1;
//cout<<"КОРНИ:"<<"z1="<<x1<<"  z2="<<x2<<"  P(z)="<<e<<"\n";getch();
   ZZ2:	;
  }//KOREN4

//**********************************************************
void KORAL(int N, double A[], complex<double> x[], int sop[]){
  // N -степень ур-я
  // A[N+1] -коэффициенты
  // x[N] -корни
  // sop[N] -прзнак комплексно-сопряженности
  //             sop=0 -действительный корень
  //             sop=1,2,... - номер пары комплексно-сопряженных корней

  double C[10], B[10];
  int M, n, i, so; 
  double u, v;
  for(n=0; n<N; n++) sop[n]=0;
  //N=N1;
  for(n=0;n<=N;n++)  B[n]=A[n]/A[N];
  if(N==1){
	  x[N-1]=-B[0]/B[1];
	  goto ZZ;
  }
  for(M=N; M>=2; M=M-2){
//     cout<<"сейчас войдем в БАРСТОУ\n";getch();
     BARSTOE(M, B, x[N-M], x[N-M+1], so);
//     cout<<"прошли  БАРСТОУ\n";getch();
     if(so==0) sop[N-M]=sop[N-M+1]=0;
     if(so==1) sop[N-M]=sop[N-M+1]=((N-M+1)+1)/2;
//cout<<M<<"\n";getch();
     u=real(-x[N-M]-x[N-M+1]);
	 v=real(x[N-M]*x[N-M+1]);

//cout<<"KORAL:M="<<M<<"\n";
//for(n=0;n<=M;n++)cout<<"B["<<n<<"]="<<B[n]<<"\n";
//cout<<"u="<<u<<"  v="<<v<<"\n";getch();

     DELEN2(M, B, u, v, n, C);
//cout<<"пройдена DELEN2\n";getch();
     for(n=0; n<=M; n++) B[n]=C[n];
//cout<<"пройдена DELEN2\n";getch();
    }//M
  if(M==1){
	  //cout<<"C[0]="<<C[0]<<"C[1]="<<C[1]<<"\n";getch();
      x[N-1]=-C[0]/C[1];}
  ZZ: ;
//cout<<"oiiu\n";getch();
    }//KORAL

//**********************************************************
void KORALkr(int N,double A[],complex<double> x[],int sop[]){
	double C[10],B[10],u,v,kr,kr1;
    int M, I, J, L, Q, j, i, n; 
    complex<double> f,y[10];

for(n=0; n<N; n++)  x[n]=complex<double>(1e10,1e10);

  for(n=N; n>=0; n--)  C[n]=A[n];

  j=0;
  for(n=1; n<N; n++)
	  if(C[0]==C[n])  j=j+1;

  if(j==N-1)  C[N-1]=C[N-1]*(1.+1e-14);
	 //при равенстве всех коэф BARSTOE не работает

  i=0;
  for(n=0;n<=N;n++){
		if(C[n]==0)  i=i+1;
		 else goto ZZ2;
  }

  M=N;

ZZ2:M=N-i;
  for(n=0;n<=M;n++)C[n]=C[n+i];
	     for(n=N-1;n>=M-1;n--)  x[n]=complex<double>(0,0);

 kr=kr1=pow(fabs(C[M]/C[0]),1./M);

 for(n=1;n<M;n++){
	 C[n]=C[n]/(kr1*C[0]);
	 kr1=kr1*kr;
 }
 C[0]=SIGN(C[0]);C[M]=1.;                          //нормировка A[0]

//cout<<"KORALkr\n"<<"M="<<M<<"\n";
//for(j=M;j>=0;j--)cout<<"C["<<j<<"]="<<C[j]<<"\n";getch();
   Q=0;
 ZZ://I=M;
 I=M-Q-1;
//cout<<"M="<<M<<"  Q="<<Q<<"  I="<<I<<"\n";
   if(I<0)goto ZZ1;
//cout<<"M="<<M<<"\n";
//for(j=M;j>=0;j--)cout<<"C["<<j<<"]="<<C[j]<<"\n";getch();
   PROIZ(I, M, C, J, B);
   if(fabs(B[0])<1e-15)goto ZZ3;
//cout<<"производная "<<I<<"\n";
//for(j=J;j>=0;j--)cout<<"B["<<j<<"]="<<B[j]<<"\n";getch();
   KORAL(J, B, y, sop);
   if(I==0){
	   for(n=0;n<J;n++)x[n]=y[n];
	   goto ZZ1;
   }
//cout<<"lkjh\n";getch();
   for(j=J-1;j>=0;j--){
	   POLY(y[j], M, C, f);
//cout<<"j="<<j<<"  POLY="<<f<<"\n";getch();
      if(sqrt(norm(f))<1e-10){
		  if(imag(y[j])==0){
			  x[M-1]=y[j];
			  DELEN1(M,C,real(y[j]),L,B);
			  M=L;
			  Q=0;
//cout<<"новые коэфф:\n";
//for(i=M;i>=0;i--)cout<<"B["<<i<<"]="<<B[i]<<"\n";
//cout<<"y["<<j<<"]="<<y[j]<<"\n";getch();
	   }//if
	if(imag(y[j])!=0){
		x[M-1]=y[j];
		x[M-2]=conj(y[j]);
//cout<<"z["<<M-1<<"]="<<x[M-1]<<"  z["<<M-2<<"]="<<x[M-2]<<"\n";getch();
	   u=real(-y[j]-conj(y[j])); 
	   v=real(y[j]*conj(y[j]));
	   DELEN2(M, C, u, v, L, B);
	   M=L;   
	   Q=0;
//cout<<"новые коэфф:\n";
//for(i=M;i>=0;i--)cout<<"B["<<i<<"]="<<B[i]<<"\n";
	   }//if
	for(n=0; n<=M; n++) C[n]=B[n];
//for(n=0;n<N;n++)cout<<"x["<<n<<"]="<<x[n]<<"\n";getch();
	goto ZZ;
	}//if
      } //j
   ZZ3:Q=Q+1;goto ZZ;
   ZZ1:
//cout<<"метка ZZ1\n";getch();
   for(n=0; n<N; n++)  x[n]=x[n]/kr;
//cout<<"КОНЕЦ  KORALkr\n\n";getch();
  }//KORALkr

//*************************************************************************
void DELEN1c(int N, complex<double> A[], complex<double> v, int&M, complex<double> B[]){ 
	int n;
   M=N-1;
   B[N-1]=A[N];
   for(n=N-2;n>=0;n--)  B[n]=A[n+1]+v*B[n+1];
   }//DELEN1c

//*************************************************************************
struct KVDR{
	int N;
	double Ko;
	complex<double> Z[8]/*корни*/;
	complex<double> a[8];
}KVDR1;

//*************************************************************************
void PRODRO(int MM, double A[], int N, double B[], KVDR&KVDR1)
 {//MM -степень знаменателя
  //A[]-коэффициенты знаменателя
  //N  -степень числителя
  //B[]-коэффициенты числителя знаменателя
  complex<double> D[110], Bc[10], Ac[10], C[10];
  int M;
  double a,b;
  int m,n,sop[10];
//cout<<"PRODRO\n"; getch();
//  KORALkr(N,B,KVDR1.Z,sop);
  KORAL(N, B, KVDR1.Z, sop);
//for(n=0;n<N;n++)cout<<"КОРЕНЬ "<<n+1<<":  z="<<KVDR1.Z[n]<<"\n";getch();
  for(n=0; n<=N; n++)  Bc[n]=complex<double>(B[n],0)/B[N];
  for(n=0; n<N; n++)  //перебор корней=столбцов
   {
	   DELEN1c(N, Bc, KVDR1.Z[n], M, C);
	   for(m=0; m<N; m++) D[n+m*(N+1)]=C[m];
  }//n
  if(MM==N){
	  KVDR1.Ko=A[N]/B[N];
  } else{
	  KVDR1.Ko=0;
  }

  for(m=0; m<N; m++){
	  if(m<=MM) a=(A[m]-B[m]*KVDR1.Ko)/B[N]; 
	  else a=0;
    D[(m+1)*(N+1)-1]=complex<double>(a,0);
//cout<<(m+1)*(N+1)-1<<"  "<<D[(m+1)*(N+1)-1]<<"\n";
    }//m

 //for(m=0;m<N;m++)
 // {for(n=0;n<N+1;n++)cout<<D[m*(N+1)+n]<<"  ";cout<<"\n";} getch();

  GAUSS1c(N,D,KVDR1.a);
  KVDR1.N=N;
cout<<"Mчис="<<MM<<"  Nзн=" <<N<<"  KVDR1.Ko="<<KVDR1.Ko<<"\n";
for(n=0; n<N; n++)cout<<"n="<<n<<"  ALF="<<KVDR1.a[n]<<"  Z="<<KVDR1.Z[n]<<"\n";
_getch();
  }//PRODRO

//**********************************************************************
//**********************************************************************
#endif __WLIB4_H
