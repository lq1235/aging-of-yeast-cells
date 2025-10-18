#include <stdio.h> 
#include <math.h>
#include <stdlib.h>
#include "stdio.h"
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define TINY 1.0e-20;

#define Dim 4
#define La 200
#define Lb 200
#define XMAX 500
#define DC 0.5
//#include <iostream>
//#include <limits> // 用于包含 std::numeric_limits 宏
double a1;
FILE *fp1, *fp2;

void mnewt(int ntrial, double x[], double tolx, double tolf);
double gasdev(long * idum);
void force(double x[], double fvec[]);
double ran1(long *idum);
double dl(double xi[], double xf[]);
double findMin(double arr[], int size);

main() 
{	
    long   r = 1, *point;
    point = &r;
    char st1[20], st2[20];
	int ntrial = 51;                 
    int i, j,  filenum, num;               
	double h_sqrt, gx[Dim], aux[Dim], xh[Dim], fvec[Dim], fxh[Dim],noise[Dim];   
	double x[Dim],  xi[Dim],xf[Dim];
	double ti, tau, tau2, var_t; 	
	double tolx = 1.0e-4, tolf = 1.0e-4;  

	double   h =0.1; 
	double  r0 = 10; 
	double  Tn =1000500000;
	double x1[Dim] = {0}, x2[Dim] = {0}; 
	double diff[Dim] = {DC, DC, DC, DC};
	double x_max[Dim] = {XMAX, XMAX, XMAX, XMAX}, x_min[Dim] = {0.0, 0.0, 0.0, 0.0};
//	double x_min[Dim]={55,0,20},x_max[Dim]={75,13,50};;  
		
	sprintf(st1, "MFPT_h=%0.3f_r=%0.2f.txt", h, r0); 
	fp1 = fopen(st1, "w+");
	//for (filenum = 1; filenum <= 8; filenum++) 
	{   
	    tau = 0;   
	   // tau2=50000; 
	//	printf("filenum:%d\n", filenum);  
	
	// a1=0.04- (filenum-1.0)*0.0001;
	// a1=0.056+ (filenum-1.0)*0.0005;
		a1=0.5;                                  
	
	       	//x1[0] = 1.586221;x1[1] = 0;x1[2] = 1.849857;x1[3] = 0;x1[4] = 1.849857;x1[5] = 0;  			
   					

	       	x2[0] = 10.15080637;x2[1] = 94.29954643;x2[2] = 272.0908633; x2[3] = 93.07199785;//  初始点 
	      	      	
	       	x1[0] =   12.21309;x1[1] =   170.60732; x1[2] = 225.94219;x1[3] = 166.11765;  //  末点 
	       	


	       	    // 	x1[0] = 30;x1[1] = 300; x1[2] = 300;x1[3] = 300;  //  末点 
		//	mnewt(ntrial,x1,tolx,tolf);    
		//	mnewt(ntrial,x2,tolx,tolf);    
		//	printf("x1=%f,%f,%f,%f,%f,%f\n", x1[0], x1[1], x1[2], x1[3], x1[4], x1[5]); 
	     //	printf("x2=%f,%f,%f,%f,%f,%f\n", x2[0], x2[1], x2[2], x2[3], x2[4], x2[5]); 
	/****
		x1[0] = 0.1483;x1[1] = 0.2027;x1[2] = 0.1891;x1[3] = 0.0721;x1[4] = 0.1902;x1[5] = 0.3268,x1[6] = 0.0748;  //有病 
		x2[0] = 1;x2[1] = 0;x2[2] = 1;x2[3] = 0, x2[4] = 0;x2[5] = 0,x2[6] = 0;         //无病 
	****/	
		printf("x1=%f,%f,%f,%f\n", x1[0], x1[1], x1[2],x1[3]); 
		printf("x2=%f,%f,%f,%f\n", x2[0], x2[1], x2[2], x2[3]); 
		
		for (i = 0; i < Dim; i++)
		{
			xi[i] = x2[i]; 
		}	
		for (i = 0; i < Dim; i++)
		{
			xf[i] = x1[i]; 	
			gx[i] = sqrt(2 * diff[i]);
		}
		sprintf(st2, "FPT_h=%0.8f_r0=%0.8f_a1=%0.5f.txt", h, r0,a1); 
		fp2 = fopen(st2, "w+"); 
	
		for (num = 1; num <= 100; num++) 
		{
			
			
			//	printf("%d\n", num); 
			for (i = 0; i < Dim; i++) 
			{
				x[i] = xi[i];        
			}
	
			h_sqrt = sqrt(h);
			for (j = 1; j <= Tn; j++) 
			{ 
				force(x, fvec);
				for (i = 0; i < Dim; i++) 
				{
					noise[i] = h_sqrt * gasdev(point);
					aux[i] = h * fvec[i] + noise[i] * gx[i];
					xh[i] = x[i] + aux[i];
				}
				force(xh, fxh);
				for (i = 0; i < Dim; i++) 
				{
					x[i] = x[i] + 0.5 * (aux[i] + h * fxh[i] + noise[i] * gx[i]);                         
					//x[i]=x[i]+fvec[i]*h+gx[i]*noise[i];//Euler method
					//xn[i]=x[i]+h*fvec[i]+sqrt(2*diff[i])*noise[i]+0.5*diffdx[i]*(pow(noise[i],2)-h);     
					if (x[i] < x_min[i])                                                                 
						x[i] = 2 * x_min[i] - x[i]; 
					else if (x[i] > x_max[i])
						x[i] = 2 * x_max[i] - x[i];
					else
						x[i] = x[i];
				}
			
				if (dl(x, xf) <= r0) 
				{                             
					ti = j * h;               
					fprintf(fp2, "%d	%0.8f\n", num, ti); 
                    printf("%d	%0.15f\n",num, ti);   
					break;                     
				}
			}
			tau = tau + ti;    
			tau2 = tau2 + ti * ti; 
			
		}
       
		tau = tau / (num -1); 
		tau2 = tau2 / (num - 1); 
		var_t = (tau2 - tau * tau) / (tau * tau); 
		printf("a1=%0.8f,tau=%0.8f,vart=%f\n", a1, tau, var_t);
		fprintf(fp1, "%0.8f	%0.8f	%f	%f	%f	%f	%f	%f	%f	%f	%f\n", a1, tau, var_t, xi[0], xi[1], xi[2], xi[3], xf[0],xf[1], xf[2], xf[3]);
	}                                                                                                                               
	fclose(fp1);
	fclose(fp2);
}

void force(double x[], double fvec[]) 
{
	
double as = 30.5;
double ah = 183;
double ah0 = 0.1;
double as0 = a1;
double beta = 3.7;//3.7或者4.6 
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

	
}
// 寻找数组中的最小值
//double findMin(double arr[], int size) {
  //  double minVal = std::numeric_limits<double>::max();
  //  for (int i = 0; i < size; ++i) {
    //    if (arr[i] < minVal) {
   //         minVal = arr[i];
   //     }
   // }
  //  return minVal;
//}



double fmax(double a, double b) 
{
	return (a > b) ? a : b;
}

double fmin(double a, double b) 
{
	return (a < b) ? a : b;
}

double dl(double xi[], double xf[]) 
{ 
	int i;
	double l = 0;
	for (i = 0; i < Dim; i++)
		l = l + (xf[i] - xi[i]) * (xf[i] - xi[i]);
	l = sqrt(l);

	return l;
}

double dot(double x[], double y[]) 
{
	int i = 0;
	double sum = 0;
	for (i = 0; i < Dim; i++)
		sum = sum + x[i] * y[i];
	return sum;
}

//Newton-Raphson methold solve the root of nonlinear equations.
void mnewt(int ntrial, double x[], double tolx, double tolf) 
{
	void lubksb(double a[Dim][Dim], int n, int indx[], double b[]);
	void ludcmp(double a[Dim][Dim], int n, int indx[], double * d);
	void fdjac(double x[], double fvec[], double df[Dim][Dim], void (*vecfunc)(double [], double []));
	int k, i, indxx[Dim] = {0};
	double errx, errf, d, a = 0;
	double pp[Dim], fvec[Dim], df[Dim][Dim] = {0};

	for (k = 0; k < ntrial; k++) 
	{
		//printf("x=%f,%f,%f,%f,%f,%f\n",x[0],x[1],x[2],x[3],x[4],x[5]);
		force(x, fvec);
		//printf("fvec1=%f,%f,%f,%f,%f,%f\n",fvec[0],fvec[1],fvec[2],fvec[3],fvec[4],fvec[5]);
		fdjac(x, fvec, df, force);
		//printf("fvec2=%f,%f,%f,%f,%f,%f\n",fvec[0],fvec[1],fvec[2],fvec[3],fvec[4],fvec[5]);
		errf = 0.0;
		for (i = 0; i < Dim; i++)
			errf += fabs(fvec[i]);
		if (errf <= tolf)
			return;
		for (i = 0; i < Dim; i++)
			pp[i] = -fvec[i];
		ludcmp(df, Dim, indxx, &d);
		lubksb(df, Dim, indxx, pp);
		errx = 0.0;
		for (i = 0; i < Dim; i++) 
		{
			errx += fabs(pp[i]);
			x[i] += pp[i];
		}
		//printf("x=%f,%f,%f,%f,%f,%f\n",x[0],x[1],x[2],x[3],x[4],x[5]);
		//getchar();
		if (errx <= tolx)  
			return;
	}
	return;
}

void ludcmp(double a[Dim][Dim], int n, int indx[], double *d) 
{
	int i, imax = 0, j, k;
	double big, dum, sum, temp;
	double vv[Dim];
	*d = 1.0;
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++)
			if ((temp = fabs(a[i][j])) > big)
				big = temp;
		//if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		if (big == 0.0) {
			big = TINY
		}//{printf("Singular matrix in routine ludcmp");exit(0);}
		vv[i] = 1.0 / big;
	}
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			if ( (dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0)
			a[j][j] = TINY;
		if (j != n - 1) {
			dum = 1.0 / (a[j][j]);
			for (i = j + 1; i < n; i++)
				a[i][j] *= dum;
		}
	}
}

void lubksb(double a[Dim][Dim], int n, int indx[], double b[]) 
{
	int i, ii = 0, ip, j;
	double sum;
	for (i = 0; i < n; i++) 
	{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii != 0)
			for (j = ii - 1; j < i; j++)
				sum -= a[i][j] * b[j];
		else if (sum != 0.0)
			ii = i + 1;
		b[i] = sum;
	}
	//printf("aa");
	for (i = n - 1; i >= 0; i--) 
	{
		sum = b[i];
		for (j = i + 1; j < n; j++)
			sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}

void fdjac(double x[], double fvec[], double df[Dim][Dim], void (*vecfunc)(double [], double [])) 
{
	int i, j;
	double h, temp, ff[Dim] = {1.0}, eps = 1.0e-6;
	for (j = 0; j < Dim; j++)
	{
		temp = x[j];
		h = eps * fabs(temp);
		if (h == 0.0)
			h = eps;
		x[j] = temp + h;
		h = x[j] - temp;
		(*vecfunc)(x, ff);
		x[j] = temp;
		for (i = 0; i < Dim; i++) 
		{
			df[i][j] = (ff[i] - fvec[i]) / h;   
		}
	}
}

double gasdev(long *idum) 
{

	static int iset = 0;
	static double gset;
	double fac, rsq, v1, v2;

	if (*idum < 0)
		iset = 0;
	if  (iset == 0) 
	{
		do 
		{
			v1 = 2.0 * ran1(idum) - 1.0;
			v2 = 2.0 * ran1(idum) - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} 
		while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq);
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	} else {
		iset = 0;
		return gset;
	}
}

double ran1(long *idum) 
{
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0 || !iy) 
	{
		if (-(*idum) < 1)
			*idum = 1;
		else *idum = -(*idum);
		for (j = NTAB + 7; j >= 0; j--) 
		{
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if (*idum < 0)
				*idum += IM;
			if (j < NTAB)
				iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k;
	if (*idum < 0)
		*idum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = *idum;
	if ((temp = AM * iy) > RNMX)
		return RNMX;
	else
		return temp;
}
