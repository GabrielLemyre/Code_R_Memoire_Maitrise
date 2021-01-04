#include <R.h>
#include <R_ext/Random.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>



void Filtering_2Regimes (int *T, int *B,double *par, double *y, double *logLikSum, double *X_t, double *X_tlag) {
 
double mu0=par[0];
double mu1=par[1];
double sigma2=par[2];

double p11=par[3];
double p22=par[4];

double p12=1-p11;
double p21=1-p22;

double p1, p2;
p1=(1-p22)/(2-p11-p22);
p2=(1-p11)/(2-p11-p22);

/*t=1*/
X_tlag[0]=p11*p1+p21*p2;
X_tlag[1]=p12*p1+p22*p2;


double *eta, *logLik;
eta = (double *) R_alloc(2*T[0], sizeof(double));
logLik = (double *) R_alloc(T[0], sizeof(double));

eta[0]=0;
eta[0]=exp(-0.5*log(2*M_PI)-0.5*log(sigma2)-0.5*pow((y[0]-mu0),2)/sigma2);
eta[1]=exp(-0.5*log(2*M_PI)-0.5*log(sigma2)-0.5*pow((y[0]-mu1),2)/sigma2);

logLik[0]=log(eta[0]*X_tlag[0]+eta[1]*X_tlag[1]);
X_t[0]=eta[0]*X_tlag[0]/exp(logLik[0]);
X_t[1]=eta[1]*X_tlag[1]/exp(logLik[0]);

for (int t=1; t < *T; t++) {

X_tlag[t*2]		=p11*X_t[(t-1)*2]+p21*X_t[(t-1)*2+1];
X_tlag[(t*2)+1]	=p12*X_t[(t-1)*2]+p22*X_t[(t-1)*2+1];

eta[t*2+0]=exp(-0.5*log(2*M_PI)-0.5*log(sigma2)-0.5*pow((y[t]-mu0),2)/sigma2);
eta[t*2+1]=exp(-0.5*log(2*M_PI)-0.5*log(sigma2)-0.5*pow((y[t]-mu1),2)/sigma2);

logLik[t]=log(eta[t*2]*X_tlag[t*2]+eta[(t*2)+1]*X_tlag[(t*2)+1]);
X_t[t*2+0]=eta[t*2+0]*X_tlag[t*2+0]/exp(logLik[t]);
X_t[t*2+1]=eta[t*2+1]*X_tlag[t*2+1]/exp(logLik[t]);

}


logLikSum[0]=0;


for (int t=*B; t < *T; t++) {
logLikSum[0]+=logLik[t];
}

}
