
function_Sim.Study<-function(i,kind,T,B,par.true){

sim<-simulation(par.true,T+B,type=kind,seed=i)
y<-sim$data.sim
pi11.true<-sim$pi11
pi22.true<-sim$pi22
S.true<-sim$S
x.sim.true<-sim$x.sim

#par(mfrow=c(2,2))
#ts.plot(y,main="Simulated data",ylab="")
#ts.plot(S.true-1,main="Simulated regimes",ylab="")
#ts.plot(pi11.true,main="Simulated p11",ylab="",ylim=c(0,1))
#ts.plot(pi22.true,main="Simulated p22",ylab="",ylim=c(0,1))
#dev.off()

####################################################################################
par<-c(par.true,0.5,0.5)
par.trasf_Const(par.trasf.inv_Const(par))
MS.est<-nlminb(start= par.trasf.inv_Const(par),filtering.single.trasf_Const, y=y,B=B,control=list(eval.max=1e6,iter.max=1e6,trace=0))

MS.par.fin<-par.trasf_Const(MS.est$par)
MS.filt<-filtering_Const(par=MS.par.fin,y,B)
MS.smooth<-smoothing_Const(MS.par.fin,X.f= attr(MS.filt,"X.t"),X.flag= attr(MS.filt,"X.tlag"))

H<-hessian(filtering.single.trasf_Const, par.trasf.inv_Const(MS.par.fin),y=y, B=B)
se<-tryCatch(stand.error.function_Const(par=MS.par.fin,hessian=H), error=function(e) return(rep(0, length(MS.par.fin))))

res.MS<-matrix(NA, ncol=3, nrow=length(MS.par.fin)+7)
res.MS[1:length(MS.par.fin),]<-cbind(MS.par.fin,se,abs(MS.par.fin/se))

res.MS[length(MS.par.fin)+1,1]<-(-MS.est$objective)
res.MS[length(MS.par.fin)+2,1]<-2*length(MS.par.fin)+2*MS.est$objective +2*length(MS.par.fin)*(length(MS.par.fin)+1)/(length(y)-B-length(MS.par.fin)-1)
res.MS[length(MS.par.fin)+3,1]<-2*MS.est$objective+length(MS.par.fin)*(log(length(y)-B)-log(2*pi))
res.MS[length(MS.par.fin)+4,1]<-MAE.function(MS.par.fin,y,attr(MS.filt,"X.tlag"),B)
res.MS[length(MS.par.fin)+5,1]<-MSE.function(MS.par.fin,y,attr(MS.filt,"X.tlag"),B)
res.MS[length(MS.par.fin)+6,1]<-MASE.function(MS.par.fin,y,attr(MS.filt,"X.tlag"),B)
res.MS[length(MS.par.fin)+7,1]<-MSSE.function(MS.par.fin,y,attr(MS.filt,"X.tlag"),B)

colnames(res.MS)<-c("Est.","SE","T test")
rownames(res.MS)<-c("mu0","mu1","sigma2","pi11","pi22","logLik","AICc","BIC","MAE","MSE","MASE","MSSE")
round(res.MS,3)

dia.res<-matrix(NA, nrow=3, ncol=7)
rownames(dia.res)<-c("MS","TVP","GAS")
colnames(dia.res)<-c("logLik","AICc","BIC","MAE","MSE","MASE","MSSE")

dia.res[1,]<-res.MS[(length(MS.par.fin)+1):(length(MS.par.fin)+7),1]


MS.smooth.series<-MS.smooth[,2]
MS.filt.p11<-rep(MS.par.fin[4],length(y))
MS.filt.p22<-rep(MS.par.fin[5],length(y))

latent.res<-matrix(NA, nrow=3, ncol=6)
rownames(latent.res)<-c("MS","TVP","GAS")
colnames(latent.res)<-c("MSE S","MAE S","MSE p11","MAE p11","MSE p22","MAE p22")
latent.res [1,1]<-mean((MS.smooth.series[(B+1):length(y)]-(S.true[(B+1):length(y)]-1))^2)
latent.res [1,2]<-mean(abs(MS.smooth.series[(B+1):length(y)]-(S.true[(B+1):length(y)]-1)))

latent.res [1,3]<-mean((MS.filt.p11[(B+1):length(y)]-pi11.true[(B+1):length(y)])^2)
latent.res [1,4]<-mean(abs(MS.filt.p11[(B+1):length(y)]-pi11.true[(B+1):length(y)]))

latent.res [1,5]<-mean((MS.filt.p22[(B+1):length(y)]-pi22.true[(B+1):length(y)])^2)
latent.res [1,6]<-mean(abs(MS.filt.p22[(B+1):length(y)]-pi22.true[(B+1):length(y)]))

round(latent.res,5)


####################################################################################
par<-c(par.true,0.5,0.5,0,0)
par.trasf_TVP(par.trasf.inv_TVP(par))
TVP.est<-nlminb(start= par.trasf.inv_TVP(par),filtering.single.trasf_TVP, y=y,B=B,control=list(eval.max=1e6,iter.max=1e6,trace=0))

TVP.par.fin<-par.trasf_TVP(TVP.est$par)
TVP.filt<-filtering_TVP(par=TVP.par.fin,y,B)
TVP.smooth<-smoothing_TVP(attr(TVP.filt,"p11"),attr(TVP.filt,"p22"),X.f= attr(TVP.filt,"X.t"),X.flag= attr(TVP.filt,"X.tlag"))

H<-hessian(filtering.single.trasf_TVP, par.trasf.inv_TVP(TVP.par.fin),y=y, B=B)
se<-tryCatch(stand.error.function_TVP(par= TVP.par.fin,hessian=H), error=function(e) return(rep(0, length(TVP.par.fin))))

res.TVP<-matrix(NA, ncol=3, nrow=length(TVP.par.fin)+7)
res.TVP[1:length(TVP.par.fin),]<-cbind(TVP.par.fin,se,abs(TVP.par.fin/se))

res.TVP[length(TVP.par.fin)+1,1]<-(-TVP.est$objective)
res.TVP[length(TVP.par.fin)+2,1]<-2*length(TVP.par.fin)+2*TVP.est$objective +2*length(TVP.par.fin)*(length(TVP.par.fin)+1)/(length(y)-B-length(TVP.par.fin)-1)
res.TVP[length(TVP.par.fin)+3,1]<-2*TVP.est$objective+length(TVP.par.fin)*(log(length(y)-B)-log(2*pi))
res.TVP[length(TVP.par.fin)+4,1]<-MAE.function(TVP.par.fin,y,attr(TVP.filt,"X.tlag"),B)
res.TVP[length(TVP.par.fin)+5,1]<-MSE.function(TVP.par.fin,y,attr(TVP.filt,"X.tlag"),B)
res.TVP[length(TVP.par.fin)+6,1]<-MASE.function(TVP.par.fin,y,attr(TVP.filt,"X.tlag"),B)
res.TVP[length(TVP.par.fin)+7,1]<-MSSE.function(TVP.par.fin,y,attr(TVP.filt,"X.tlag"),B)


colnames(res.TVP)<-c("Est.","SE","T test")
rownames(res.TVP)<-c("mu0","mu1","sigma2","pi11","pi22","A1","A2","logLik","AICc","BIC","MAE","MSE","MASE","MSSE")
round(res.TVP,3)

dia.res[2,]<-res.TVP[(length(TVP.par.fin)+1):(length(TVP.par.fin)+7),1]

TVP.smooth.series<-TVP.smooth[,2]
TVP.filt.p11<-attr(TVP.filt,"p11")
TVP.filt.p22<-attr(TVP.filt,"p22")

latent.res [2,1]<-mean((TVP.smooth.series[(B+1):length(y)]-(S.true[(B+1):length(y)]-1))^2)
latent.res [2,2]<-mean(abs(TVP.smooth.series[(B+1):length(y)]-(S.true[(B+1):length(y)]-1)))

latent.res [2,3]<-mean((TVP.filt.p11[(B+1):length(y)]-pi11.true[(B+1):length(y)])^2)
latent.res [2,4]<-mean(abs(TVP.filt.p11[(B+1):length(y)]-pi11.true[(B+1):length(y)]))

latent.res [2,5]<-mean((TVP.filt.p22[(B+1):length(y)]-pi22.true[(B+1):length(y)])^2)
latent.res [2,6]<-mean(abs(TVP.filt.p22[(B+1):length(y)]-pi22.true[(B+1):length(y)]))


round(dia.res,3)
round(latent.res,3)


###################################
par<-c(par.true,0.5,0.5,0,0,0.9,0.9)

requireNamespace("numDeriv")
requireNamespace("statmod")

par.trasf_GAS(par.trasf.inv_GAS(par))
GAS.est<-nlminb(start= par.trasf.inv_GAS(par),filtering.single.trasf_GAS, y=y,B=B,control=list(eval.max=1e6,iter.max=1e6,trace=0))

GAS.par.fin<-par.trasf_GAS(GAS.est$par)
GAS.filt<-filtering_GAS(par=GAS.par.fin,y,B)
GAS.smooth<-smoothing_GAS(attr(GAS.filt,"p11"),attr(GAS.filt,"p22"),X.f= attr(GAS.filt,"X.t"),X.flag= attr(GAS.filt,"X.tlag"))


H<-hessian(filtering.single.trasf_GAS, par.trasf.inv_GAS(GAS.par.fin),y=y, B=B)
se<-tryCatch(stand.error.function_GAS(par= GAS.par.fin,hessian=H), error=function(e) return(rep(0, length(GAS.par.fin))))



res.GAS<-matrix(NA, ncol=3, nrow=length(GAS.par.fin)+7)
res.GAS[1:length(GAS.par.fin),]<-cbind(GAS.par.fin,se,abs(GAS.par.fin/se))

res.GAS[length(GAS.par.fin)+1,1]<-(-GAS.est$objective)
res.GAS[length(GAS.par.fin)+2,1]<-2*length(GAS.par.fin)+2*GAS.est$objective +2*length(GAS.par.fin)*(length(GAS.par.fin)+1)/(length(y)-B-length(GAS.par.fin)-1)
res.GAS[length(GAS.par.fin)+3,1]<-2*GAS.est$objective+length(GAS.par.fin)*(log(length(y)-B)-log(2*pi))
res.GAS[length(GAS.par.fin)+4,1]<-MAE.function(GAS.par.fin,y,attr(GAS.filt,"X.tlag"),B)
res.GAS[length(GAS.par.fin)+5,1]<-MSE.function(GAS.par.fin,y,attr(GAS.filt,"X.tlag"),B)
res.GAS[length(GAS.par.fin)+6,1]<-MASE.function(GAS.par.fin,y,attr(GAS.filt,"X.tlag"),B)
res.GAS[length(GAS.par.fin)+7,1]<-MSSE.function(GAS.par.fin,y,attr(GAS.filt,"X.tlag"),B)


colnames(res.GAS)<-c("Est.","SE","T test")
rownames(res.GAS)<-c("mu0","mu1","sigma2","pi11","pi22","A1","A2","B1","B2","logLik","AICc","BIC","MAE","MSE","MASE","MSSE")
round(res.GAS,3)

dia.res[3,]<-res.GAS[(length(GAS.par.fin)+1):(length(GAS.par.fin)+7),1]

GAS.smooth.series<-GAS.smooth[,2]
GAS.filt.p11<-attr(GAS.filt,"p11")
GAS.filt.p22<-attr(GAS.filt,"p22")

latent.res [3,1]<-mean((GAS.smooth.series[(B+1):length(y)]-(S.true[(B+1):length(y)]-1))^2)
latent.res [3,2]<-mean(abs(GAS.smooth.series[(B+1):length(y)]-(S.true[(B+1):length(y)]-1)))

latent.res [3,3]<-mean((GAS.filt.p11[(B+1):length(y)]-pi11.true[(B+1):length(y)])^2)
latent.res [3,4]<-mean(abs(GAS.filt.p11[(B+1):length(y)]-pi11.true[(B+1):length(y)]))

latent.res [3,5]<-mean((GAS.filt.p22[(B+1):length(y)]-pi22.true[(B+1):length(y)])^2)
latent.res [3,6]<-mean(abs(GAS.filt.p22[(B+1):length(y)]-pi22.true[(B+1):length(y)]))


round(dia.res,3)
round(latent.res,3)



###################################################################################
if(kind=="TVP-AR"){

par<-c(par.true,0.5,0.5,0,0)
par.trasf_TVPXExo(par.trasf.inv_TVPXExo(par))
TVP.XExo.est<-nlminb(start= par.trasf.inv_TVPXExo(par),filtering.single.trasf_TVPXExo, y=y,B=B, X_Exo = x.sim.true,control=list(eval.max=1e6,iter.max=1e6,trace=0))

TVP.XExo.par.fin<-par.trasf_TVPXExo(TVP.XExo.est$par)
TVP.XExo.filt<-filtering_TVPXExo(par=TVP.XExo.par.fin,y,B,x.sim.true)
TVP.XExo.smooth<-smoothing_TVPXExo(attr(TVP.XExo.filt,"p11"),attr(TVP.XExo.filt,"p22"),X.f= attr(TVP.XExo.filt,"X.t"),X.flag= attr(TVP.XExo.filt,"X.tlag"))

H<-hessian(filtering.single.trasf_TVPXExo, par.trasf.inv_TVPXExo(TVP.XExo.par.fin),y=y, B=B, X_Exo = x.sim.true)
se<-tryCatch(stand.error.function_TVPXExo(par=TVP.XExo.par.fin,hessian=H), error=function(e) return(rep(NA, length(TVP.XExo.par.fin))))

res.TVP.XExo<-matrix(NA, ncol=3, nrow=length(TVP.XExo.par.fin)+7)
res.TVP.XExo[1:length(TVP.XExo.par.fin),]<-cbind(TVP.XExo.par.fin,se,abs(TVP.XExo.par.fin/se))

res.TVP.XExo[length(TVP.XExo.par.fin)+1,1]<-(-TVP.XExo.est$objective)
res.TVP.XExo[length(TVP.XExo.par.fin)+2,1]<--2*length(TVP.XExo.par.fin)+2*TVP.XExo.est$objective +2*length(TVP.XExo.par.fin)*(length(TVP.XExo.par.fin)+1)/(length(y)-B-length(TVP.XExo.par.fin)-1)
res.TVP.XExo[length(TVP.XExo.par.fin)+3,1]<-2*TVP.XExo.est$objective+length(TVP.XExo.par.fin)*(log(length(y)-B)-log(2*pi))
res.TVP.XExo[length(TVP.XExo.par.fin)+4,1]<-MAE.function(TVP.XExo.par.fin,y,attr(TVP.XExo.filt,"X.tlag"),B)
res.TVP.XExo[length(TVP.XExo.par.fin)+5,1]<-MSE.function(TVP.XExo.par.fin,y,attr(TVP.XExo.filt,"X.tlag"),B)
res.TVP.XExo[length(TVP.XExo.par.fin)+6,1]<-MASE.function(TVP.XExo.par.fin,y,attr(TVP.XExo.filt,"X.tlag"),B)
res.TVP.XExo[length(TVP.XExo.par.fin)+7,1]<-MSSE.function(TVP.XExo.par.fin,y,attr(TVP.XExo.filt,"X.tlag"),B)


colnames(res.TVP.XExo)<-c("Est.","SE","T test")
rownames(res.TVP.XExo)<-c("mu0","mu1","sigma2","pi11","pi22","A1","A2","logLik","AICc","BIC","MAE","MSE","MASE","MSSE")

round(res.TVP.XExo,3)

dia.res=rbind(dia.res,"TVP-X"=NA)
dia.res[4,]<-res.TVP.XExo[(length(TVP.XExo.par.fin)+1):(length(TVP.XExo.par.fin)+7),1]


TVP.XExo.smooth.series<-TVP.XExo.smooth[,2]
TVP.XExo.filt.p11<-attr(TVP.XExo.filt,"p11")
TVP.XExo.filt.p22<-attr(TVP.XExo.filt,"p22")

latent.res=rbind(latent.res,"TVP-X"=NA)

latent.res [4,1]<-mean((TVP.XExo.smooth.series[(B+1):length(y)]-(S.true[(B+1):length(y)]-1))^2)
latent.res [4,2]<-mean(abs(TVP.XExo.smooth.series[(B+1):length(y)]-(S.true[(B+1):length(y)]-1)))

latent.res [4,3]<-mean((TVP.XExo.filt.p11[(B+1):length(y)]-pi11.true[(B+1):length(y)])^2)
latent.res [4,4]<-mean(abs(TVP.XExo.filt.p11[(B+1):length(y)]-pi11.true[(B+1):length(y)]))

latent.res [4,5]<-mean((TVP.XExo.filt.p22[(B+1):length(y)]-pi22.true[(B+1):length(y)])^2)
latent.res [4,6]<-mean(abs(TVP.XExo.filt.p22[(B+1):length(y)]-pi22.true[(B+1):length(y)]))


round(dia.res,3)
round(latent.res,3)

}



names.dia<-expand.grid(rownames(dia.res),colnames(dia.res))
vec.dia.res<-c(dia.res)
names(vec.dia.res)<-paste0(as.character(names.dia$Var1),"_",as.character(names.dia$Var2))

names.latent<-expand.grid(rownames(latent.res),colnames(latent.res))
vec.latent.res<-c(latent.res)
names(vec.latent.res)<-paste0(as.character(names.latent$Var1),"_",as.character(names.latent$Var2))

res<-c(vec.dia.res, vec.latent.res)
names(res)<-names(c(vec.dia.res, vec.latent.res))
return(res)
}

