
filtering_Const<-function(par, y,B)
{
y<-as.numeric(y)
T<-length(y)

temp<-.C("Filtering_2Regimes","T"=as.integer(T),"B"=as.integer(B),"par"=as.numeric(par), "y"=as.double(y),"logLikSum"=double(1),"X_t"=double(2*T),"X_tlag"=double(2*T),NAOK=TRUE,PACKAGE="HMMGAS")
res<-(-temp$logLikSum)

attr(res,"X.t")<-matrix(temp$X_t,ncol=2,byrow=T)
attr(res,"X.tlag")<-matrix(temp$X_tlag,ncol=2,byrow=T)

res
}


smoothing_Const<-function(par,X.f, X.flag)
{

p11<-par[4]
p12<-1-p11

p22<-par[5]
p21<-1-p22

P<-matrix(c(p11,p12,p21,p22),2,2,byrow=T)

T<-nrow(X.f)
X<-matrix(0, ncol=2, nrow=T)

t=T
X[t,]<-X.f[t,]

t=T-1
for (t in ((T-1):1)) X[t,]<-(P%*%(X[t+1,]/X.flag[t+1,]))*X.f[t,]

return("X.s"=X)
}



par.trasf_Const<-function(par)
{
par.t<-par
par.t[3]<-exp(par[3])
par.t[c(4,5)]<-exp(par[c(4,5)])/(1+exp(par[c(4,5)]))
return(par.t)}

par.trasf.inv_Const<-function(par)
{par.t<-par
par.t[3]<-log(par[3])
par.t[c(4,5)]<-log(par[c(4,5)])-log(1-par[c(4,5)])
return(par.t)}


stand.error.function_Const<-function(par,hessian)
{
H<-hessian

G.d<-rep(1, length(par))
G.d[c(3)]<-1/par[c(3)]
G.d[c(4,5)]<-1/(par[c(4,5)])+1/(1-par[c(4,5)])
G<-matrix(0,length(par),length(par))
diag(G)<-G.d

se<-diag(solve(G)%*%solve(H)%*%solve(t(G)))^0.5

se	
}


filtering.single.trasf_Const<-function(par,y,B){
par.t<-par.trasf_Const(par)
l<-filtering_Const(par.t,y,B)
l[1:1]}







