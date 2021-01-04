
MAE.function<-function(par,data,X.t,B)
{
mu.a<-par[1]
mu.b<-par[2]
sigma2<-par[3]

T<-length(data)
SFE<-rep(0,T)
for (t in B:T) SFE[t]<-(data[t]-X.t[t,1]*mu.a-X.t[t,2]*mu.b)
mean(abs(SFE[B:T]))	}

MSE.function<-function(par,data,X.t,B)
{
mu.a<-par[1]
mu.b<-par[2]
sigma2<-par[3]

T<-length(data)
SFE<-rep(0,T)
for (t in B:T) SFE[t]<-(data[t]-X.t[t,1]*mu.a-X.t[t,2]*mu.b)
mean((SFE[B:T])^2)	}


MASE.function<-function(par,data,X.t,B)
{
mu.a<-par[1]
mu.b<-par[2]
sigma2<-par[3]

T<-length(data)
SFE<-rep(0,T)
for (t in B:T) {
mean.t<-X.t[t,1]*mu.a+X.t[t,2]*mu.b
var.t<-sigma2+X.t[t,1]*(mu.a-mean.t)^2+X.t[t,2]*(mu.b-mean.t)^2
SFE[t]<-(data[t]-mean.t)/sqrt(var.t)}
mean(abs(SFE[B:T]))	}

MSSE.function<-function(par,data,X.t,B)
{
mu.a<-par[1]
mu.b<-par[2]
sigma2<-par[3]

T<-length(data)
SFE<-rep(0,T)
for (t in B:T) {
mean.t<-X.t[t,1]*mu.a+X.t[t,2]*mu.b
var.t<-sigma2+X.t[t,1]*(mu.a-mean.t)^2+X.t[t,2]*(mu.b-mean.t)^2
SFE[t]<-(data[t]-mean.t)/sqrt(var.t)}
mean((SFE[B:T])^2)	}

