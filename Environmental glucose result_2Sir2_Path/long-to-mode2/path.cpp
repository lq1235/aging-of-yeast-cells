#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M 41                                                 
#define Dim 2
#define Min_num 2
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define DC 0.05                       

#include <iostream>
#include <limits> 


FILE *fp;
double diff[Dim][Dim],invd[Dim][Dim],invdx[Dim][Dim][Dim];
double e_eff,lambda=100.0;
double vx[Min_num][Dim]={   0.82595,0.971274,      	1.71797 ,	0.0144378 }; 

double saddle[Dim]={0,0}, diag[Dim][Dim]={1,0,0,1};
double x_min[Dim]={0,0},x_max[Dim]={2.5,2.5};
double a1=1.1;

double veff(double x[]);
void path_mc(double xi[],double xf[],double path[M][Dim]);
double veff_min(double xi[]);
void diffx(double x[]);
void force(double xp[],double fvec[]);
void fdjac(double x[], double fvec[], double df[Dim][Dim],void (*vecfunc)(double [], double []));
double dl(double inv[Dim][Dim],double xi[],double xf[]);
double dot(double x[],double y[]);
double dshj(double xi[],double xf[]);
double dt(double xi[], double xf[]);
double xmax(double a,double b);
double xmin(double a,double b);
double time(double x[M][Dim]);
double shj_lambda(double x[M][Dim]);
double shj(double x[M][Dim]);

void path_mc(double xi[],double xf[],double x_path[M][Dim]);
double veff_min(double x[]);
double findMin(double arr[], int size);
main()
{
    int n,ndiff,fixpoint,i,j,k,nk=1,num,grid[Dim];
    double  itgl[Min_num],transition_time,temp=0,pb=0,phi,l=0,t=10;
    double c=0,e_old,amplitude=1,ds=0,deltat=0,smin=0,e1=0,e2=0,x1[Dim],x2[Dim];
    double ps[Min_num][M][Dim],psr[Min_num][M][Dim],x_path[M][Dim]={0};

    char st[20];
    num=21;
    
    
    e_eff = -xmin(veff(vx[0]),veff(vx[1]));  
    //e_eff=0.0;
   
	printf("e1=%f,e2=%f,e_eff=%f\n",veff(vx[0]),veff(vx[1]),e_eff);
	
 
    
	    fixpoint=1;
		//define the staring and ending position
		
		for(i=0;i<Dim;i++)
		{
			ps[fixpoint-1][0][i]=vx[0][i];            
			ps[fixpoint-1][M-1][i]=vx[1][i];         
		}
                path_mc(ps[fixpoint-1][0],ps[fixpoint-1][M-1],x_path);

		for(j=1;j<M-1;j++)
			for(i=0;i<Dim;i++)
				ps[fixpoint-1][j][i]=x_path[j][i];
				//printf("hehe\n");
				itgl[fixpoint-1]=shj(ps[fixpoint-1]);
  
	smin=itgl[0];
            printf("phi=%f,fixpoint=%d\n",itgl[fixpoint-1],fixpoint);
            sprintf(st,"path.txt");           
            //sprintf(st,"path%0.2f_a1_%0.1f.txt",DC,a1);           
	    printf("num=%d\n",num);
            if((fp=fopen(st,"w"))==NULL)
            {
                printf("Cannot open file. \n");
                exit(0);
            }
         
            for(i=0;i<M;i++)
            {

			fprintf(fp,"%f	%f	%f	%f\n",1.0*i,ps[fixpoint-1][i][0],ps[fixpoint-1][i][1],ps[fixpoint-1][i][2]);  

	        }

            fclose(fp);
            if((fp=fopen("shj.txt","a+"))==NULL)  
            {
                printf("Cannot open file. \n");
                exit(0);
            }
		transition_time=time(ps[fixpoint-1]);
		fprintf(fp,"num	M	a1	lambda	e_eff	x0	x1	x2	Shj	time\n");
		fprintf(fp,"%d	%d	%0.1f	%0.1e	%0.2f	%0.2f	%0.2f	%0.2f	%f	%f\n",num,M,a1,lambda,e_eff,ps[fixpoint-1][M-1][0],ps[fixpoint-1][M-1][1],ps[fixpoint-1][M-1][2],itgl[fixpoint-1],transition_time);
        fclose(fp);
}
double veff(double x[])
{
    int i,j,k;
    void force(double x[],double fvec[]);
	void fdjac(double x[], double fvec[], double df[Dim][Dim],void (*vecfunc)(double [], double []));
	double m=0.01,bb=0,cc=0,vf;
	double fvec[Dim],df[Dim][Dim]={0};
	diffx(x);
	force(x,fvec);
	fdjac(x,fvec,df,force);
    for(i=0;i<Dim;i++)
		for(j=0;j<Dim;j++)
		{
			bb=bb+invd[i][j]*fvec[i]*fvec[j];
			for(k=0;k<Dim;k++)
				cc=cc+diff[j][k]*invd[i][j]*df[i][k]+diff[j][k]*fvec[i]*invdx[k][i][j];
		}
   
     vf=0.25*bb+0.5*cc;  
	  
	return vf;
} 

void path_mc(double xi[],double xf[],double path[M][Dim])
{
	int ii,i,j,k,n,m,pj,mk,acceptance_num=0,cycle=1,times=5000;
	int inttau,tau,pos_start,pos_end,pos_ex,n_vari=3,Steps=1000,si;//+++++
	long ra=1,*point;
	double itgl=0,itglr=0,l=0,r,t=1.0e-2,t_init,tfactor1=0.9,rho=0.05,tfactor2=0.7;
	double S,S1,delta,delta_x,avg_old,avg_new,fluc_old,fluc_new,dE;
	double psr[M][Dim],amplitude={0},ran1(long *idum),xc[Dim];
	point=&ra;
	for(i=0;i<Dim;i++){
		psr[0][i]=path[0][i]=xi[i];
		psr[M-1][i]=path[M-1][i]=xf[i];
	}
	
	

	for(i=1;i<=M-2;i++){
		for(n=0;n<Dim;n++){
			path[i][n]=psr[i][n]=psr[0][n]+(psr[M-1][n]-psr[0][n])*i/(M-1.0);
			
			
		}
	}

    itgl=shj(path);
    printf("shj_straight=%f\n",itgl);
        printf("i	amplitude	t	acceptance_ratio	shj	penalty\n");
		t_init=t;
        amplitude=0.0001; 
	for(i=1;i<=Steps;i++)
	{
		acceptance_num=0;
		for(si=1;si<=times;si++)
		{
			S = 0;
			S1 = 0;
			do
			{
				pos_start=rand()%M;
				pos_end=rand()%M;
				if (pos_start > pos_end) {
				pos_ex = pos_end;
				pos_end = pos_start;
				pos_start = pos_ex;
				}
			} while(pos_end-pos_start<2);

				for (tau=pos_start+1; tau<=pos_end; tau++)
				{
					S = S + dshj(path[tau-1],path[tau]);
					if (tau<pos_end)
					{
						for(n=0;n<Dim;n++)
						{
							inttau = rand();
							delta = amplitude*rand()/(double)RAND_MAX;
							delta_x = delta*(x_max[n]-x_min[n])/(double)M;
							//delta_x = delta/(double)M;
							if (inttau%n_vari == 0) psr[tau][n] = path[tau][n];
							if (inttau%n_vari == 1 && path[tau][n]+delta_x <= x_max[n]) psr[tau][n] = path[tau][n] + delta_x;
							if (inttau%n_vari == 2 && path[tau][n]-delta_x >= x_min[n]) psr[tau][n] = path[tau][n] - delta_x;
						}
					}
					for(m=0;m<Dim;m++) xc[m]=(psr[tau-1][m]+psr[tau][m])/2.0;
					if(e_eff+veff(xc)>=0) {
						S1 = S1 + dshj(psr[tau-1],psr[tau]);
					}
					else {
						S1 = S1 + 10000.0;
					}
				}
				avg_old = 0;
				avg_new = 0;
				for (k=pos_start+1; k<=pos_end; k++)
				{
					avg_old = avg_old + dl(diag,path[k-1],path[k]);
					avg_new = avg_new + dl(diag,psr[k-1],psr[k]);
				}
				avg_old = avg_old/(pos_end-pos_start);
				avg_new = avg_new/(pos_end-pos_start);
				fluc_old=0;
				fluc_new=0;
				for (k=pos_start+1; k<=pos_end; k++)
				{
					fluc_old = fluc_old + pow(dl(diag,path[k-1],path[k])-avg_old,2);
					fluc_new = fluc_new + pow(dl(diag,psr[k-1],psr[k])-avg_new,2);
				}
				S = S + lambda*fluc_old;
				S1 = S1 + lambda*fluc_new;
				if (S>S1)
				{
					for(k=pos_start+1;k<pos_end;k++)
					{
						for (n=0;n<Dim;n++) path[k][n] = psr[k][n];
					}
					S = S1;
					acceptance_num++;
				}
				for(k=pos_start+1;k<pos_end;k++)
					{
						for (n=0;n<Dim;n++) psr[k][n] = path[k][n];
					}
		}
		itgl=shj(path);
		rho=acceptance_num*1.0/times;
		printf("%d	%f	%0.2e	%0.3f	%f	%f\n",i,amplitude,t,rho,itgl,shj_lambda(path)-itgl);
		t=t*tfactor1;
	}
}

double veff_min(double xi[])
{
	int ii,i,j,k,n,m,pj,mk,acceptance_num=0,times=10000;
	int inttau,n_vari=2,Steps=100,si;
	double v0,v,l=0,r,t=1.0e-5,t_init,tfactor1=0.9,rho=0.05,tfactor2=0.7;
	double S,S1,delta,delta_x,x[Dim],x_new[Dim];
	double psr[M][Dim],amplitude={0},ran1(long *idum),xc[Dim];
	
	for (n=0;n<Dim;n++) x[n] = x_new[n] = xi[n];
    v0=veff(xi);
	t_init=t;
    amplitude=0.001;
	for(i=1;i<=Steps;i++)
	{
		acceptance_num=0;
		for(si=0;si<=times;si++)
		{
			S = veff(x);
			for(n=0;n<Dim;n++)
				{
					inttau = rand();
					delta = amplitude*rand()/(double)RAND_MAX;
					delta_x = delta*(x_max[n]-x_min[n]);
					if (inttau%n_vari == 0 && x[n]+delta_x <= x_max[n]) x_new[n] = x[n] + delta_x;
					if (inttau%n_vari == 1 && x[n]-delta_x >= x_min[n]) x_new[n] = x[n] - delta_x;
				}
			S1 = veff(x_new);
			if (S>S1)
			{
				for (n=0;n<Dim;n++) x[n] = x_new[n];
				S = S1;
				acceptance_num++;
			}
		}
		v=veff(x);
		rho=acceptance_num*1.0/times;
		t=t*tfactor1; 	
	}
	for(j=0;j<Dim;j++) printf("x[%d]=%f	",j,x[j]);
	for(j=0;j<Dim;j++) xi[j] = x[j];
	printf("\n");
	printf("v0=%f	vf=%f\n",v0,v);
	return v;
}

void diffx(double x[])
{
	int i,j,k;
	
	diff[0][0]=DC;
	diff[0][1]=0;
//	diff[0][2]=0;


	
	diff[1][0]=0;
    diff[1][1]=DC;
   // diff[1][2]=0;


    
  //  diff[2][0]=0;
  //  diff[2][1]=0;
  //  diff[2][2]=DC;

    

    
    // inverse of diffusion coefficient
	invd[0][0]=1.0/diff[0][0];
	invd[0][1]=0;
//invd[0][2]=0;

	
	
	invd[1][0]=0;
	invd[1][1]=1.0/diff[1][1];
//	invd[1][2]=0;

	
//	invd[2][0]=0;
//	invd[2][1]=0;
//	invd[2][2]=1.0/diff[2][2];




    // d D ^-1 / dx 
    for(i=0;i++;i<=Dim){
    	for(j=0;j++;j<=Dim){
    		for(k=0;k++;k<=Dim){
    			invdx[i][j][k]=0;
    		}
    	}
    }

}

void force(double x[],double fvec[])
{

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

}


double findMin(double arr[], int size) {
    double minVal = std::numeric_limits<double>::max();
    for (int i = 0; i < size; ++i) {
        if (arr[i] < minVal) {
            minVal = arr[i];
        }
    }
    return minVal;
}


void fdjac(double x[], double fvec[], double df[Dim][Dim],void (*vecfunc)(double [], double []))
{
	int i,j;
	double h,temp,ff[Dim]={1.0},eps=1.0e-6;

	for (j=0;j<Dim;j++) {
		temp=x[j];
		h=eps*fabs(temp);
		if (h == 0.0) h=eps;
		x[j]=temp+h;
		h=x[j]-temp;
		(*vecfunc)(x,ff);
		x[j]=temp;
		for (i=0;i<Dim;i++){
			df[i][j]=(ff[i]-fvec[i])/h;
		}
	}
}

double dl(double inv[Dim][Dim],double xi[],double xf[])
{
    int i,j;
    double l=0;
    for(i=0;i<Dim;i++)
		for(j=0;j<Dim;j++)
			l=l+inv[i][j]*(xf[i]-xi[i])*(xf[j]-xi[j]);
    l=sqrt(l);
    return l;
}

double dot(double x[],double y[])
{
    int i=0;
    double sum=0;
    for(i=0;i<Dim;i++)
        sum=sum+x[i]*y[i];
    return sum;
}

double dshj(double xi[],double xf[])
{
	int i,j;
    double xc[Dim],ds=0,vi,fdl=0,dx[Dim],fd[Dim];
    double check_rest(double x[]),veff(double xi[]),dot(double x[],double y[]);
    double fvec[Dim];
	for(i=0;i<Dim;i++) xc[i]=(xf[i]+xi[i])/2;
	force(xc,fvec);
	vi=veff(xc);
	for(i=0;i<Dim;i++) dx[i]=xf[i]-xi[i];
	for(i=0;i<Dim;i++)
		for(j=0;j<Dim;j++)
			fdl=fdl+invd[i][j]*fvec[j]*dx[i];
//	ds=dl(invd,xi,xf)*sqrt(vi+e_eff)-fdl/2.0;
		//ds=dl(invd,xi,xf)*sqrt(vi+e_eff);//-fdl/2.0;
		ds=dl(invd,xi,xf)*sqrt(vi+e_eff)-fdl/2.0;  
	
    return ds;
}

double dt(double xi[], double xf[])
{
	int i;
	double xc[Dim],vi,dt=0;
	double veff(double xi[]);
	for(i=0;i<Dim;i++) xc[i]=(xf[i]+xi[i])/2;
	//for(i=0;i<Dim;i++) xc[i]=xi[i];
	vi=veff(xc);
	dt=dl(diag,xi,xf)/sqrt(4*(vi+e_eff));
	return dt;
}

double xmax(double a,double b)
{
	return (a>b)?a:b;
}

double xmin(double a,double b)
{
	return (a<b)?a:b;
}

double time(double x[M][Dim])
{
	int i;
	double sum_t=0;
	double dt(double xi[],double xf[]);
	for(i=0;i<=M-2;i++)
		sum_t=sum_t+dt(x[i],x[i+1]);
	return sum_t;
}

double shj_lambda(double x[M][Dim])
{
	int i;
	double s=0,lave=0;
	double dl(double inv[Dim][Dim],double xi[],double xf[]),dshj(double xi[],double xf[]);
	lave=0;	
	for(i=0;i<=M-2;i++)
	{
		diffx(x[i]);
		lave=lave+dl(diag,x[i],x[i+1]);
	}
	lave=lave/(M-1);
	for(i=0;i<=M-2;i++)
	{
		s=s+dshj(x[i],x[i+1])+lambda*pow(dl(diag,x[i],x[i+1])-lave,2);
	}
	return s;
}

double shj(double x[M][Dim])
{
	int i;
	double s=0,lave=0;
	double dl(double inv[Dim][Dim],double xi[],double xf[]),dshj(double xi[],double xf[]);

	for(i=0;i<=M-2;i++)
	{
		s=s+dshj(x[i],x[i+1]);
	}
	return s;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
#undef M
#undef Dim
#undef Gn
#undef Min_num

