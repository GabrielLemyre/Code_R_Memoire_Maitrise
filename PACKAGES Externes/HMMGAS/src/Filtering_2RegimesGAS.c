#include <R.h>
#include <R_ext/Random.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>


void Filtering_2RegimesGAS (int *T, int *B,double *par, double *y, double *logLikSum, double *X_t, double *X_tlag, double *p11, double *p22, double *nodes, double *weights, double *mu_weights, double *sigma_weights) {

double mu0=par[0];
double mu1=par[1];
double sigma0=par[2];
double sigma1=par[3];

double omega1_LR=log(par[4]/(1-par[4]));
double omega2_LR=log(par[5]/(1-par[5]));

double A1=par[6];
double A2=par[7];

double B1=par[8];
double B2=par[9];

double omega1=omega1_LR;
double omega2=omega2_LR;


double *p1, *p2,*f1,*f2,*score_SCAL, *I, *eta,*logLik;

f1= (double*) R_alloc( T[0] , sizeof(double) ) ;
f2= (double*) R_alloc( T[0] , sizeof(double) ) ;
p1= (double*) R_alloc( T[0] , sizeof(double) ) ;
p2= (double*) R_alloc( T[0] , sizeof(double) ) ;
eta = (double *) R_alloc(2*T[0], sizeof(double));
logLik = (double *) R_alloc(T[0], sizeof(double));
score_SCAL = (double *) R_alloc(2*T[0], sizeof(double));
I = (double *) R_alloc(T[0], sizeof(double));



f1[0]=omega1_LR;
f2[0]=omega2_LR;


p11[0]=1/(1+exp(-f1[0]));
p22[0]=1/(1+exp(-f2[0]));

// p1[0]=(1-p22[0])/(2-p11[0]-p22[0]);
// p2[0]=(1-p11[0])/(2-p11[0]-p22[0]);


X_tlag[0]=p11[0]*p1[0]+(1-p22[0])*p2[0];
X_tlag[1]=(1-p11[0])*p1[0]+p22[0]*p2[0];

eta[0]=exp(-0.5*log(2*M_PI)-0.5*log(sigma0)-0.5*pow((y[0]-mu0),2)/sigma0);
eta[1]=exp(-0.5*log(2*M_PI)-0.5*log(sigma1)-0.5*pow((y[0]-mu1),2)/sigma1);

logLik[0]=log(eta[0]*X_tlag[0]+eta[1]*X_tlag[1]);
X_t[0]=eta[0]*X_tlag[0]/exp(logLik[0]);
X_t[1]=eta[1]*X_tlag[1]/exp(logLik[0]);

double d1, d2, den;
double I_star;
I_star=0;

I[0]=0;

for (int j=0; j < 30; j++) {


d1=exp(-0.5*log(2*M_PI)-0.5*log(sigma0)-0.5*pow((nodes[j]-mu0),2)/sigma0);
d2=exp(-0.5*log(2*M_PI)-0.5*log(sigma1)-0.5*pow((nodes[j]-mu1),2)/sigma1);

den=d1*(p11[0]*p1[0]+(1-p22[0])*p2[0])+d2*((1-p22[0])*p1[0]+p22[0]*p2[0]);

if(den ==0){
I_star=0;
} else{
// printf("den at 1 is %f \n", den);

I_star=(pow((d1-d2),2)/den)*pow(2*M_PI*sigma_weights[0]*sigma_weights[0], 0.5)*exp(pow(nodes[j]-mu_weights[0],2)/(2*sigma_weights[0]*sigma_weights[0]));
}

I[0]=I[0]+I_star*weights[j];
}



double S, *g;
g= (double*) R_alloc( T[0]*2 , sizeof(double) ) ;


S=(eta[0]-eta[1])/exp(logLik[0]);
g[0]=	 p1[0]*p11[0]*(1-p11[0]);
g[1]=	-p2[0]*p22[0]*(1-p22[0]);

double g_mod;
g_mod=pow((g[0]*g[0]+g[1]*g[1]), 0.5);

if(g_mod==0) {
g[0]=0;
g[1]=0;
} else{

g[0]=g[0]/g_mod;
g[1]=g[1]/g_mod;

}

double S_star;

if(I[0]==0){
S_star=0;
} else{
S_star=S/pow(I[0],0.5);
}

score_SCAL[0]= S_star*g[0];
score_SCAL[1]= S_star*g[1];

f1[1]=omega1+ A1*score_SCAL[0]+ B1*(f1[0]-omega1);
f2[1]=omega2+ A2*score_SCAL[1]+ B2*(f2[0]-omega2);


p11[1]=1e-10+(1-2*1e-10)/(1+exp(-f1[1]));
p22[1]=1e-10+(1-2*1e-10)/(1+exp(-f2[1]));



for (int t=1; t < *T; t++) {


X_tlag[t*2]		=p11[t]		*X_t[(t-1)*2]+	(1-p22[t])	*X_t[(t-1)*2+1];
X_tlag[(t*2)+1]	=(1-p11[t])	*X_t[(t-1)*2]+	p22[t]		*X_t[(t-1)*2+1];

eta[t*2+0]=exp(-0.5*log(2*M_PI)-0.5*log(sigma0)-0.5*pow((y[t]-mu0),2)/sigma0);
eta[t*2+1]=exp(-0.5*log(2*M_PI)-0.5*log(sigma1)-0.5*pow((y[t]-mu1),2)/sigma1);

logLik[t]=log(eta[t*2]*X_tlag[t*2]+eta[(t*2)+1]*X_tlag[(t*2)+1]);
X_t[t*2+0]=eta[t*2+0]*X_tlag[t*2+0]/exp(logLik[t]);
X_t[t*2+1]=eta[t*2+1]*X_tlag[t*2+1]/exp(logLik[t]);

logLik[t]=log(eta[t*2]*X_tlag[t*2]+eta[t*2+1]*X_tlag[t*2+1]);
X_t[t*2+1]=eta[t*2+1]*X_tlag[t*2+1]/exp(logLik[t]);
X_t[t*2]=1-X_t[t*2+1];


I_star=0;
I[t]=0;

for (int j=0; j < 30; j++) {

d1=exp(-0.5*log(2*M_PI)-0.5*log(sigma0)-0.5*pow((nodes[j]-mu0),2)/sigma0);
d2=exp(-0.5*log(2*M_PI)-0.5*log(sigma1)-0.5*pow((nodes[j]-mu1),2)/sigma1);

den=d1*(p11[t]*X_t[(t-1)*2]+(1-p22[t])*X_t[(t-1)*2+1])+d2*((1-p22[t])*X_t[(t-1)*2]+p22[t]*X_t[(t-1)*2+1]);

if(den ==0){
I_star=0;
} else{
// printf("den at %d is %f \n",t, den);

I_star=(pow((d1-d2),2)/den)*pow(2*M_PI*sigma_weights[0]*sigma_weights[0], 0.5)*exp(pow(nodes[j]-mu_weights[0],2)/(2*sigma_weights[0]*sigma_weights[0]));
}


I[t]=I[t]+I_star*weights[j];
}


S=(eta[t*2]-eta[t*2+1])/exp(logLik[t]);
g[t*2]=	 	X_t[(t-1)*2]*p11[t]*(1-p11[t]);
g[t*2+1]=	-X_t[(t-1)*2+1]*p22[t]*(1-p22[t]);

double g_mod;
g_mod=pow((g[t*2]*g[t*2]+g[t*2+1]*g[t*2+1]), 0.5);

if(g_mod==0) {
g[t*2]=0;
g[t*2+1]=0;
} else{

g[t*2]=g[t*2]/g_mod;
g[t*2+1]=g[t*2+1]/g_mod;

}

double S_star;

if(I[t]==0){
S_star=0;
} else{
S_star=S/pow(I[t],0.5);
}

score_SCAL[t*2]= S_star*g[t*2];
score_SCAL[t*2+1]= S_star*g[t*2+1];

if((t<(T[0]-1))){

f1[t+1]=omega1+ A1*score_SCAL[t*2]		+B1*(f1[t]-omega1);
f2[t+1]=omega2+ A2*score_SCAL[t*2+1]	+B2*(f2[t]-omega2);

p11[t+1]=1e-10+(1-2*1e-10)/(1+exp(-f1[t+1]));
p22[t+1]=1e-10+(1-2*1e-10)/(1+exp(-f2[t+1]));
}


}



logLikSum[0]=0;

for (int t=*B; t < *T; t++) {
logLikSum[0]+=logLik[t];

}

}
