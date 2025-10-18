   /************ trajectory and Probability************/
#include <math.h>
#include <stdlib.h>
#include "stdio.h"

#define TINY 1.0e-20;
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define Dim 4
#define La 200
#define Lb 200
#define XMAX 500
#define DC 0.5

FILE *fp1,*fp2;
char st[20],st1[20],st2[20];
double a1; 

main()
{
	int i,j,m,k,index1,index2,A,B,filenum;
	int grid[Dim],p_num[La*Lb]={0},p[La][Lb]={0};
	double s=0,eps=1.0e-4,f_cyc[La][Lb]={0},f_cdk1[La][Lb]={0};
	double fvec[Dim]={0},x[Dim]={24.2676387075027,	161.244264460529,	324.816463610110,	206.962295436469},xc[2],xi[2],xf[2];
	void rk2(double h,double x[],double fvec[],double xmin[],double xmax[],long *point);
	double x_min[Dim],x_max[Dim],tau,tau0;
	double h=0.01,x_old[Dim];
	
	double iter=1.0,Tn=2.0e10,Tni=2.0e6,Tnf,n_write=1.0e4;
	long   r=1,*point;
	
	Tnf = Tn+Tni;
	point=&r;
	for (i=0;i<Dim;i++) {
		x_max[i]=XMAX;
		x_min[i]=0.0;
	//	x[i]=1.5;
		}
	//  Òª Þ¸ÄµÄ²   	
		char a_1[20]="D1"; 
		
		
for(filenum=1;filenum<=21;filenum++)
{
	a1=3.3+ (filenum-1.0)*0.1; 
//	a1=0.1; 
	
//	a1[0][1]=a;
	for(A=0;A<La;A++)
	{
		for(B=0;B<Lb;B++)
		{
			p[A][B]= 0.0;
			f_cyc[A][B] = 0.0;
			f_cdk1[A][B] = 0.0;
		}
	}
	//a1 = 0.5 + (filenum-1.0)*0.5;
	index1 = 3;
	index2 = 4;
	printf("filenum:%d\n",filenum);
	//sprintf(st1,"bistable_a%0.1f.txt",a);
	sprintf(st1,"bistable_%0.5f.txt",a1);
	if((fp1=fopen(st1,"w+"))==NULL){
		printf("Cannot open file. \n");exit(0);}
	tau = 0;
	tau0 = 0;
	j = 1;
	for(iter=1.0;iter<=Tnf;iter=iter+1)
	{ 
//		if(iter==Tni+1)
//		{
//			fprintf(fp1,"%0.0f	%f",iter-Tni-1,tau0);
//			//printf("%e %f",iter-Tni-1,tau0);
//			for (i=0;i<Dim;i++){
//				fprintf(fp1,"	%f",x[i]);
//				//printf("	%f",x[i]);
//			}
//			fprintf(fp1,"\n");
//			//printf("\n");	
//		}
		rk2(h,x,fvec,x_min,x_max,point);	
		if(iter>=Tni&iter<=Tni+n_write)
		{
			fprintf(fp1,"%0.0f	%f",iter-Tni,tau0);
			tau0=tau0+h;
			//printf("%e %f",iter-Tni,tau0);
			for (i=0;i<Dim;i++){
				fprintf(fp1,"	%f",x[i]);
				//printf("	%f",x[i]);
			}
			fprintf(fp1,"\n");
			//printf("\n");
			//getchar();
		}
		if(j%10000000==0){
			printf("iter=%e\n",iter);
			j=0;
		}
		j++;
		xc[0]=x[index1-1]; xi[0]=x_min[index1-1]; xf[0]=x_max[index1-1];
		xc[1]=x[index2-1]; xi[1]=x_min[index2-1]; xf[1]=x_max[index2-1];
		if(iter>=Tni) 
		{
			A = (int)((xc[0]-xi[0])*La/(xf[0]-xi[0]));
			B = (int)((xc[1]-xi[1])*Lb/(xf[1]-xi[1]));
			p[A][B] = p[A][B] + 1;
			f_cyc[A][B] = f_cyc[A][B] + fvec[index1-1];
			f_cdk1[A][B] = f_cdk1[A][B] + fvec[index2-1];
		}
		tau=tau+h;
	}
		//sprintf(st2,"p1_a%0.1f.txt",a);
		sprintf(st2,"pp_%s=%0.5f.txt",a_1,a1);
		if ((fp2=fopen(st2,"w+"))==NULL){
			printf("Cannot open file. \n");exit(0);}
		for(A=0;A<La;A++){
			for(B=0;B<La;B++){
					if(p[A][B]==0) fprintf(fp2,"%d	%d	0	0.0	0.0\n",A,B);
					else fprintf(fp2,"%d	%d	%d	%e	%e\n",A,B,p[A][B],f_cyc[A][B]/p[A][B],f_cdk1[A][B]/p[A][B]);
			}
		}
		fclose(fp1);
	fclose(fp2);
}	

}

void force(double x[],double fvec[])
{
	
	
double as = 30.5;
double ah = 183;
double ah0 = 0.1;
double as0 = a1;
double beta = 3.7;//3.7»òÕß4.6 
double dm = 0.3;
double dh = 3.8;
double ds = 0.2;
double Kh = 326;
double Ks = 185;
double n1 = 3;
double n2 = 4.8;


/**
mS=x(1);
mH=x(2);
S=x(3);
H=x(4);
**/

 fvec[0] = as0+as*pow(x[3], n1)/(pow(Kh, n1)+pow(x[3], n1))-dm*x[0];

  fvec[1] = ah0+ah*pow(Ks, n2)/(pow(Ks, n2)+pow(x[2], n2))-dm*x[1];

    fvec[2] = beta*x[0]-ds*x[2];

      fvec[3] = beta*x[1]-dh*x[3];


	
/**	
double beta1 = 0.01, beta2 = 0.01;
//double n1 = 3, n2 = 3;
double K_S = 0.52, K_H = 0.62;
double d_S = 0.1, d_H = 0.3;
double gamma1 = 1;

double a_St = 0.0515, K_St = 0.03, k_St = 0.02036, b_St = 1.6207;
double a_Ht = 0.126, K_Ht = 0.04, k_Ht = 0.06, b_Ht = 1.999;

double total_S,total_H;
double D1=a1, gamma2 = 0.6;

total_S=1.2*(a_St *(pow(K_St, 3)/ (pow(K_St, 3) + pow(D1, 3) ))- k_St * D1 + b_St);
total_H=(a_Ht * pow(K_Ht, 3) / (pow(K_Ht, 3) + pow(D1, 3))) - k_Ht * D1 + b_Ht;


	fvec[0]=(beta1 + pow(x[0], 3) / (pow(K_S, 3) + pow(x[0], 3))) * (total_S - x[0]) - (d_S + gamma1 * x[1]) * x[0];
	
	
	fvec[1]=(beta2 + pow(x[1], 3) / (pow(K_H, 3) + pow(x[1], 3))) * (total_H - x[1]) - (d_H + gamma2 * x[0]) * x[1];
	
	**/
	
	
	
	
	
	
}

void rk2(double h,double x[],double fvec[],double xmin[],double xmax[],long *point)
{
	int i,j,k,indx[Dim];
	//long   r=1,*point;
	double xh[Dim],fxh[Dim],sqrt_h;
	double diff[Dim],gx[Dim],k1[Dim],xn[Dim],noise[Dim];
	void force(double x[],double fvec[]);
	double gasdev(long *idum);
	//point=&r;
	for (i=0;i<Dim;i++) {
		diff[i]=DC;
		gx[i]=sqrt(2*diff[i]);
		}	
	//while(1)
	//{
		sqrt_h = sqrt(h);
		force(x,fvec);
		for(i=0;i<Dim;i++)
		{
			noise[i]=sqrt_h*gx[i]*gasdev(point);
			k1[i]=h*fvec[i]+noise[i];
			xh[i]=x[i]+k1[i];
		}
		force(xh,fxh);
		for(i=0;i<Dim;i++){
			xn[i]=x[i]+0.5*(k1[i]+h*fxh[i]+noise[i]);// Heun method
			if(xn[i]<xmin[i]) x[i]=2*xmin[i]-xn[i];//reflecting boundary condition
			else if(xn[i]>xmax[i]) x[i]=2*xmax[i]-xn[i]; 
			else x[i]=xn[i];
		}	
	//}
}

double gasdev(long *idum)
{
	double ran1(long *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if (*idum < 0) iset=0;
	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}

double ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
