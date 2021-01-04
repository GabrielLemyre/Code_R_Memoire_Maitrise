
filtering_TVP<-function(par, y,B)
{
y<-as.numeric(y)
T<-length(y)



temp<-.C("Filtering_2RegimesTVP","T"=as.integer(T),"B"=as.integer(B),"par"=as.numeric(par), "y"=as.double(y),"logLikSum"=double(1),"X_t"=double(2*T),"X_tlag"=double(2*T),"p11"=double(T),"p22"=double(T),NAOK=TRUE,PACKAGE="HMMGAS")
res<-(-temp$logLikSum)

attr(res,"X.t")<-matrix(temp$X_t,ncol=2,byrow=T)
attr(res,"X.tlag")<-matrix(temp$X_tlag,ncol=2,byrow=T)
attr(res,"p11")<-temp$p11
attr(res,"p22")<-temp$p22

res
}


smoothing_TVP<-function(p11,p22,X.f, X.flag)
{
p12<-1-p11
p21<-1-p22

T<-nrow(X.f)
X<-matrix(0, ncol=2, nrow=T)

t=T
X[t,]<-X.f[t,]

t=T-1
for (t in ((T-1):1)) 
{
P<-matrix(c(p11[t+1],p12[t+1],p21[t+1],p22[t+1]), ncol=2, byrow=T)
X[t,]<-(P%*%(X[t+1,]/X.flag[t+1,]))*X.f[t,]}

return("X.s"=X)}

par.trasf_TVP<-function(par)
{
par.t<-par
par.t[3]<-exp(par[3])
par.t[c(4,5)]<-exp(par[c(4,5)])/(1+exp(par[c(4,5)]))
return(par.t)}

par.trasf.inv_TVP<-function(par)
{par.t<-par
par.t[3]<-log(par[3])
par.t[c(4,5)]<-log(par[c(4,5)])-log(1-par[c(4,5)])
return(par.t)}


stand.error.function_TVP<-function(par,hessian)
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


filtering.single.trasf_TVP<-function(par,y,B){
par.t<-par.trasf_TVP(par)
l<-filtering_TVP(par.t,y,B)
l[1:1]}



