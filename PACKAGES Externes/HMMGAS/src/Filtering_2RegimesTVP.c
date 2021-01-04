#include <R.h>
#include <R_ext/Random.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>


void Filtering_2RegimesTVP (int *T, int *B,double *par, double *y, double *logLikSum, double *X_t, double *X_tlag, double *p11, double *p22) {
 
double mu0=par[0];
double mu1=par[1];
double sigma2=par[2];

double omega1_LR=log(par[3]/(1-par[3])); 
double omega2_LR=log(par[4]/(1-par[4]));

double A1=par[5];
double A2=par[6];

double omega1=omega1_LR*(1-A1);
double omega2=omega2_LR*(1-A2);


double *eta, *logLik;
eta = (double *) R_alloc(2*T[0], sizeof(double));
logLik = (double *) R_alloc(T[0], sizeof(double));


double *p1, *p2,*f1,*f2;

f1= (double*) R_alloc( T[0] , sizeof(double) ) ;
f2= (double*) R_alloc( T[0] , sizeof(double) ) ;
p1= (double*) R_alloc( T[0] , sizeof(double) ) ;
p2= (double*) R_alloc( T[0] , sizeof(double) ) ;


f1[0]=omega1_LR;
f2[0]=omega2_LR;


p11[0]=1/(1+exp(-f1[0]));
p22[0]=1/(1+exp(-f2[0]));

p1[0]=(1-p22[0])/(2-p11[0]-p22[0]);
p2[0]=(1-p11[0])/(2-p11[0]-p22[0]);


X_tlag[0]=p11[0]*p1[0]+(1-p22[0])*p2[0];
X_tlag[1]=(1-p11[0])*p1[0]+p22[0]*p2[0];

eta[0]=exp(-0.5*log(2*M_PI)-0.5*log(sigma2)-0.5*pow((y[0]-mu0),2)/sigma2);
eta[1]=exp(-0.5*log(2*M_PI)-0.5*log(sigma2)-0.5*pow((y[0]-mu1),2)/sigma2);

logLik[0]=log(eta[0]*X_tlag[0]+eta[1]*X_tlag[1]);
X_t[0]=eta[0]*X_tlag[0]/exp(logLik[0]);
X_t[1]=eta[1]*X_tlag[1]/exp(logLik[0]);

f1[1]=omega1;
f2[1]=omega2;


p11[1]=1/(1+exp(-f1[1]));
p22[1]=1/(1+exp(-f2[1]));


for (int t=1; t < *T; t++) {


X_tlag[t*2]		=p11[t]		*X_t[(t-1)*2]+	(1-p22[t])	*X_t[(t-1)*2+1];
X_tlag[(t*2)+1]	=(1-p11[t])	*X_t[(t-1)*2]+	p22[t]		*X_t[(t-1)*2+1];

eta[t*2+0]=exp(-0.5*log(2*M_PI)-0.5*log(sigma2)-0.5*pow((y[t]-mu0),2)/sigma2);
eta[t*2+1]=exp(-0.5*log(2*M_PI)-0.5*log(sigma2)-0.5*pow((y[t]-mu1),2)/sigma2);

logLik[t]=log(eta[t*2]*X_tlag[t*2]+eta[(t*2)+1]*X_tlag[(t*2)+1]);
X_t[t*2+0]=eta[t*2+0]*X_tlag[t*2+0]/exp(logLik[t]);
X_t[t*2+1]=eta[t*2+1]*X_tlag[t*2+1]/exp(logLik[t]);

if((t<(T[0]-1))){

f1[t+1]=omega1+ A1*y[t];
f2[t+1]=omega2+ A2*y[t];

p11[t+1]=1/(1+exp(-f1[t+1]));
p22[t+1]=1/(1+exp(-f2[t+1]));

}



}



logLikSum[0]=0;


for (int t=*B; t < *T; t++) {
logLikSum[0]+=logLik[t];

}

}
