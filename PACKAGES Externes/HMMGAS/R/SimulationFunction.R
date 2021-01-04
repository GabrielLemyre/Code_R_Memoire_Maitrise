sim.ar<-function(par)
{
mu1<-par[1]
sigma2<-par[2]
mean.y<-mu1
rnorm(1, mean.y,sqrt(sigma2))}


simulation<-function(par, size,type,seed=1)
{

set.seed(seed)
mu1<-par[1]
mu2<-par[2]
sigma2<-par[3]

par1<-c(mu1,sigma2)
par2<-c(mu2,sigma2)
par0<-c(0,0)
T<-size
x.sim=NULL

if(all(type!=c("TVP","TVP-AR"))){

if(type=="Constant"){
	pi1<-rep(0.95,T)
	pi2<-rep(0.85,T)
	}
	
	

if(type=="Break2"){
	pi1<-pi2<-rep(0,T)
	for(t in 1:T) {
		pi1[t]<-0.8-0.6*(t<(T/2))#Break
		pi2[t]<-0.2+0.6*(t<(T/2))} #Break
		}
	

if(type=="Break3"){
	pi1<-pi2<-rep(0,T)
	for(t in 1:T) {
		pi1[t]<-0.6-0.2*(t<(T/2))#Break
		pi2[t]<-0.4+0.2*(t<(T/2))} #Break
		}
	



if(type=="SlowSine"){
	pi1<-pi2<-rep(0,T)
	for(t in 1:T) {
		pi1[t]<-0.5+0.45*cos(2*pi*t/(T/2)) #SlowSine
		pi2[t]<-0.5-0.45*cos(2*pi*t/(T/2))}
		}


if(type=="Sine"){
	pi1<-pi2<-rep(0,T)
	for(t in 1:T) {
		pi1[t]<-0.5+0.45*cos(2*pi*t/(T/4)) #Sine
		pi2[t]<-0.5-0.45*cos(2*pi*t/(T/4))
		}
		}


if(type=="FastSine"){
	pi1<-pi2<-rep(0,T)
	for(t in 1:T)  pi1[t]<-0.5+0.45*cos(2*pi*t/(T/10)) #FastSine
	for(t in 1:T)  pi2[t]<-0.5-0.45*cos(2*pi*t/(T/10)) #FastSine
	}



##########
S<-data.sim<-rep(0,T)
t=1
p1<-(1-pi2[t])/(2-pi2[t]-pi1[t]) #P[S_0=1]
p2<-(1-pi1[t])/(2-pi2[t]-pi1[t]) #P[S_0=2]




S[t]<-sample(c(1,2),1,prob=c(p1,p2))
if(S[t]==1) data.sim[t]<-sim.ar(par1)
if(S[t]==2) data.sim[t]<-sim.ar(par2)



for (t in 2:T){

if(S[t-1]==1) S[t]<-sample(c(1,2),1,prob=c(pi1[t],(1-pi1[t])))
if(S[t-1]==2) S[t]<-sample(c(1,2),1,prob=c((1-pi2[t]),pi2[t]))

if((S[t]==1)) data.sim[t]<-sim.ar(par1)
if((S[t]==2)) data.sim[t]<-sim.ar(par2)

}
}else{

if(all(type==c("TVP"))){
#TVP CASE

pi1<-pi2<-rep(0,T)
S<-data.sim<-rep(0,T)

t=1

pi1[t]<-0.5 #TVP
pi2[t]<-0.5 #TVP

p1<-(1-pi2[t])/(2-pi2[t]-pi1[t]) #P[S_0=1]
p2<-(1-pi1[t])/(2-pi2[t]-pi1[t]) #P[S_0=2]

S[t]<-sample(c(1,2),1,prob=c(p1,p2))
if(S[t]==1) data.sim[t]<-sim.ar(par1)
if(S[t]==2) data.sim[t]<-sim.ar(par2)

for(t in 2:T){
	pi1[t]<-1/(1+exp(-1.2*data.sim[t-1]))#TVP
	pi2[t]<-1/(1+exp(0.3*data.sim[t-1])) #TVP
	
	p1<-(1-pi2[t])/(2-pi2[t]-pi1[t]) #P[S_0=1]
	p2<-(1-pi1[t])/(2-pi2[t]-pi1[t]) #P[S_0=2]
	
	S[t]<-sample(c(1,2),1,prob=c(p1,p2))
	
	if(S[t]==1) data.sim[t]<-sim.ar(par1)
	if(S[t]==2) data.sim[t]<-sim.ar(par2)
}	
}else{
	

#TVP-AR CASE
pi1<-pi2<-rep(0,T)
S<-data.sim<-x.sim<-rep(0,T)

t=1
pi1[t]<-0.5 #TVP-AR
pi2[t]<-0.5 #TVP-AR

x.sim[t]=rnorm(1)


p1<-(1-pi2[t])/(2-pi2[t]-pi1[t]) #P[S_0=1]
p2<-(1-pi1[t])/(2-pi2[t]-pi1[t]) #P[S_0=2]

S[t]<-sample(c(1,2),1,prob=c(p1,p2))
if(S[t]==1) data.sim[t]<-sim.ar(par1)
if(S[t]==2) data.sim[t]<-sim.ar(par2)

for(t in 2:T){
	x.sim[t]=0.99*x.sim[t-1]+sqrt(1-0.99*0.99)*rnorm(1)
	pi1[t]<-1/(1+exp(-0.7*x.sim[t]))#TVP-AR
	pi2[t]<-1/(1+exp(0.7*x.sim[t])) #TVP-AR
	
	p1<-(1-pi2[t])/(2-pi2[t]-pi1[t]) #P[S_0=1]
	p2<-(1-pi1[t])/(2-pi2[t]-pi1[t]) #P[S_0=2]
	
	S[t]<-sample(c(1,2),1,prob=c(p1,p2))
	
	if(S[t]==1) data.sim[t]<-sim.ar(par1)
	if(S[t]==2) data.sim[t]<-sim.ar(par2)
}	


	
}	
}


return(list("data.sim"=data.sim,"S"=S,"pi11"=pi1,"pi22"=pi2,"x.sim"=x.sim))}


