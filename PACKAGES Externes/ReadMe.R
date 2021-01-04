
# ##Install HMMGAS package

# # Create a new folder where to copy/paste HMMGAS_2.0.tar.gz
# # On MacOS X, open Terminal and locate the folder where you saved HMMGAS files

# # Run the following command to install the package
# # R CMD install  -l lib HMMGAS_2.0.tar.gz
# #To install into the library tree lib, use R CMD INSTALL -l lib pkgs. 
# #This prepends lib to the library path for duration of the install, so required packages in the installation #directory will be found (and used in preference to those in other libraries).

# Open R

rm(list=ls())
library(HMMGAS)

##Simulation##
par.true<-c(-1,1,0.5) #c(mu1, mu2, sigma2)

type="Sine"

# type = "TVP-AR"
# type = "Constant"
# type = "Break2"
# type = "Break3"
# type = "SlowSine"
# type = "Sine"
# type = "FastSine"
# type = "TVP"
    
size<-1000
B=200 #Burn-in
seed=1
    
sim <- simulation(par.true, size+B, type = type, seed = seed)

y <- sim$data.sim
pi11.true <- sim$pi11
pi22.true <- sim$pi22
S.true <- sim$S
x.sim.true <- sim$x.sim

##Estimation: Time constant model##
par <- c(par.true, 0.5, 0.5)

MS.est <- nlminb(start = par.trasf.inv_Const(par), filtering.single.trasf_Const, 
                   y = y, B = B, control = list(eval.max = 1e+06, iter.max = 1e+06, trace = 0))
                   
MS.est$convergence #0 means that we reached convergence
 
MS.par.fin <- par.trasf_Const(MS.est$par)

MS.par.fin

#Compare estimate parameters vs true
round(cbind(par.true ,MS.par.fin[1:3]),4)

H <- hessian(filtering.single.trasf_Const, par.trasf.inv_Const(MS.par.fin), y = y, B = B)
se <- tryCatch(stand.error.function_Const(par = MS.par.fin, 
        hessian = H), error = function(e) return(rep(0, length(MS.par.fin))))



MS.filt <- filtering_Const(par = MS.par.fin, y, B)

MS.smooth <- smoothing_Const(MS.par.fin, X.f = attr(MS.filt, "X.t"), X.flag = attr(MS.filt, "X.tlag"))

    res.MS <- matrix(NA, ncol = 3, nrow = length(MS.par.fin) + 
        7)
    res.MS[1:length(MS.par.fin), ] <- cbind(MS.par.fin, se, abs(MS.par.fin/se))
    res.MS[length(MS.par.fin) + 1, 1] <- (-MS.est$objective)
    res.MS[length(MS.par.fin) + 2, 1] <- 2 * length(MS.par.fin) + 
        2 * MS.est$objective + 2 * length(MS.par.fin) * (length(MS.par.fin) + 
        1)/(length(y) - B - length(MS.par.fin) - 1)
    res.MS[length(MS.par.fin) + 3, 1] <- 2 * MS.est$objective + 
        length(MS.par.fin) * (log(length(y) - B) - log(2 * pi))
    res.MS[length(MS.par.fin) + 4, 1] <- MAE.function(MS.par.fin, 
        y, attr(MS.filt, "X.tlag"), B)
    res.MS[length(MS.par.fin) + 5, 1] <- MSE.function(MS.par.fin, 
        y, attr(MS.filt, "X.tlag"), B)
    res.MS[length(MS.par.fin) + 6, 1] <- MASE.function(MS.par.fin, 
        y, attr(MS.filt, "X.tlag"), B)
    res.MS[length(MS.par.fin) + 7, 1] <- MSSE.function(MS.par.fin, 
        y, attr(MS.filt, "X.tlag"), B)
    colnames(res.MS) <- c("Est.", "SE", "T test")
    rownames(res.MS) <- c("mu0", "mu1", "sigma2", "pi11", "pi22", 
        "logLik", "AICc", "BIC", "MAE", "MSE", "MASE", "MSSE")
    round(res.MS, 3)
  



MS.smooth.series <- MS.smooth[, 2]
MS.filt.p11 <- rep(MS.par.fin[4], length(y))
MS.filt.p22 <- rep(MS.par.fin[5], length(y))



##Estimation: TVP Model##
par <- c(par.true, 0.5, 0.5, 0, 0)
TVP.est <- nlminb(start = par.trasf.inv_TVP(par), filtering.single.trasf_TVP, 
        y = y, B = B, control = list(eval.max = 1e+06, iter.max = 1e+06, 
            trace = 0))

    TVP.par.fin <- par.trasf_TVP(TVP.est$par)
    TVP.filt <- filtering_TVP(par = TVP.par.fin, y, B)
    TVP.smooth <- smoothing_TVP(attr(TVP.filt, "p11"), attr(TVP.filt, 
        "p22"), X.f = attr(TVP.filt, "X.t"), X.flag = attr(TVP.filt, 
        "X.tlag"))
    H <- hessian(filtering.single.trasf_TVP, par.trasf.inv_TVP(TVP.par.fin), 
        y = y, B = B)
    se <- tryCatch(stand.error.function_TVP(par = TVP.par.fin, 
        hessian = H), error = function(e) return(rep(0, length(TVP.par.fin))))
    res.TVP <- matrix(NA, ncol = 3, nrow = length(TVP.par.fin) + 
        7)
    res.TVP[1:length(TVP.par.fin), ] <- cbind(TVP.par.fin, se, 
        abs(TVP.par.fin/se))
    res.TVP[length(TVP.par.fin) + 1, 1] <- (-TVP.est$objective)
    res.TVP[length(TVP.par.fin) + 2, 1] <- 2 * length(TVP.par.fin) + 
        2 * TVP.est$objective + 2 * length(TVP.par.fin) * (length(TVP.par.fin) + 
        1)/(length(y) - B - length(TVP.par.fin) - 1)
    res.TVP[length(TVP.par.fin) + 3, 1] <- 2 * TVP.est$objective + 
        length(TVP.par.fin) * (log(length(y) - B) - log(2 * pi))
    res.TVP[length(TVP.par.fin) + 4, 1] <- MAE.function(TVP.par.fin, 
        y, attr(TVP.filt, "X.tlag"), B)
    res.TVP[length(TVP.par.fin) + 5, 1] <- MSE.function(TVP.par.fin, 
        y, attr(TVP.filt, "X.tlag"), B)
    res.TVP[length(TVP.par.fin) + 6, 1] <- MASE.function(TVP.par.fin, 
        y, attr(TVP.filt, "X.tlag"), B)
    res.TVP[length(TVP.par.fin) + 7, 1] <- MSSE.function(TVP.par.fin, 
        y, attr(TVP.filt, "X.tlag"), B)
    colnames(res.TVP) <- c("Est.", "SE", "T test")
    rownames(res.TVP) <- c("mu0", "mu1", "sigma2", "pi11", "pi22", 
        "A1", "A2", "logLik", "AICc", "BIC", "MAE", "MSE", "MASE", 
        "MSSE")
    round(res.TVP, 3)
    
    
    
    
    TVP.smooth.series <- TVP.smooth[, 2]
    TVP.filt.p11 <- attr(TVP.filt, "p11")
    TVP.filt.p22 <- attr(TVP.filt, "p22")
    
    
    
    
 ##Estimation: GAS Model##
   
    par <- c(par.true, 0.5, 0.5, 0, 0, 0.9, 0.9)
    library(MASS)
    library(statmod)
    
    par.trasf_GAS(par.trasf.inv_GAS(par))
    
    GAS.est <- nlminb(start = par.trasf.inv_GAS(par), filtering.single.trasf_GAS, 
        y = y, B = B, control = list(eval.max = 1e+06, iter.max = 1e+06, 
            trace = 0))
    
    GAS.par.fin <- par.trasf_GAS(GAS.est$par)
    
    GAS.filt <- filtering_GAS(par = GAS.par.fin, y, B)
    
    GAS.smooth <- smoothing_GAS(attr(GAS.filt, "p11"), attr(GAS.filt, 
        "p22"), X.f = attr(GAS.filt, "X.t"), X.flag = attr(GAS.filt, 
        "X.tlag"))
    
    H <- hessian(filtering.single.trasf_GAS, par.trasf.inv_GAS(GAS.par.fin), 
        y = y, B = B)
    
    se <- tryCatch(stand.error.function_GAS(par = GAS.par.fin, 
        hessian = H), error = function(e) return(rep(0, length(GAS.par.fin))))
    
    
    res.GAS <- matrix(NA, ncol = 3, nrow = length(GAS.par.fin) + 
        7)
    res.GAS[1:length(GAS.par.fin), ] <- cbind(GAS.par.fin, se, 
        abs(GAS.par.fin/se))
    res.GAS[length(GAS.par.fin) + 1, 1] <- (-GAS.est$objective)
    res.GAS[length(GAS.par.fin) + 2, 1] <- 2 * length(GAS.par.fin) + 
        2 * GAS.est$objective + 2 * length(GAS.par.fin) * (length(GAS.par.fin) + 
        1)/(length(y) - B - length(GAS.par.fin) - 1)
    res.GAS[length(GAS.par.fin) + 3, 1] <- 2 * GAS.est$objective + 
        length(GAS.par.fin) * (log(length(y) - B) - log(2 * pi))
    res.GAS[length(GAS.par.fin) + 4, 1] <- MAE.function(GAS.par.fin, 
        y, attr(GAS.filt, "X.tlag"), B)
    res.GAS[length(GAS.par.fin) + 5, 1] <- MSE.function(GAS.par.fin, 
        y, attr(GAS.filt, "X.tlag"), B)
    res.GAS[length(GAS.par.fin) + 6, 1] <- MASE.function(GAS.par.fin, 
        y, attr(GAS.filt, "X.tlag"), B)
    res.GAS[length(GAS.par.fin) + 7, 1] <- MSSE.function(GAS.par.fin, 
        y, attr(GAS.filt, "X.tlag"), B)
    colnames(res.GAS) <- c("Est.", "SE", "T test")
    rownames(res.GAS) <- c("mu0", "mu1", "sigma2", "pi11", "pi22", 
        "A1", "A2", "B1", "B2", "logLik", "AICc", "BIC", "MAE", 
        "MSE", "MASE", "MSSE")
    round(res.GAS, 3)
   
   
   
   
   
    GAS.smooth.series <- GAS.smooth[, 2]
    GAS.filt.p11 <- attr(GAS.filt, "p11")
    GAS.filt.p22 <- attr(GAS.filt, "p22")
   
   
   
    ##Estimation: TVP with autoregressive process##

   
    if (type == "TVP-AR") {
        par <- c(par.true, 0.5, 0.5, 0, 0)
        par.trasf_TVPXExo(par.trasf.inv_TVPXExo(par))
        TVP.XExo.est <- nlminb(start = par.trasf.inv_TVPXExo(par), 
            filtering.single.trasf_TVPXExo, y = y, B = B, X_Exo = x.sim.true, 
            control = list(eval.max = 1e+06, iter.max = 1e+06, 
                trace = 0))
        TVP.XExo.par.fin <- par.trasf_TVPXExo(TVP.XExo.est$par)
        TVP.XExo.filt <- filtering_TVPXExo(par = TVP.XExo.par.fin, 
            y, B, x.sim.true)
        TVP.XExo.smooth <- smoothing_TVPXExo(attr(TVP.XExo.filt, 
            "p11"), attr(TVP.XExo.filt, "p22"), X.f = attr(TVP.XExo.filt, 
            "X.t"), X.flag = attr(TVP.XExo.filt, "X.tlag"))
        H <- hessian(filtering.single.trasf_TVPXExo, par.trasf.inv_TVPXExo(TVP.XExo.par.fin), 
            y = y, B = B, X_Exo = x.sim.true)
        se <- tryCatch(stand.error.function_TVPXExo(par = TVP.XExo.par.fin, 
            hessian = H), error = function(e) return(rep(NA, 
            length(TVP.XExo.par.fin))))
        res.TVP.XExo <- matrix(NA, ncol = 3, nrow = length(TVP.XExo.par.fin) + 
            7)
        res.TVP.XExo[1:length(TVP.XExo.par.fin), ] <- cbind(TVP.XExo.par.fin, 
            se, abs(TVP.XExo.par.fin/se))
        res.TVP.XExo[length(TVP.XExo.par.fin) + 1, 1] <- (-TVP.XExo.est$objective)
        res.TVP.XExo[length(TVP.XExo.par.fin) + 2, 1] <- -2 * 
            length(TVP.XExo.par.fin) + 2 * TVP.XExo.est$objective + 
            2 * length(TVP.XExo.par.fin) * (length(TVP.XExo.par.fin) + 
                1)/(length(y) - B - length(TVP.XExo.par.fin) - 
                1)
        res.TVP.XExo[length(TVP.XExo.par.fin) + 3, 1] <- 2 * 
            TVP.XExo.est$objective + length(TVP.XExo.par.fin) * 
            (log(length(y) - B) - log(2 * pi))
        res.TVP.XExo[length(TVP.XExo.par.fin) + 4, 1] <- MAE.function(TVP.XExo.par.fin, 
            y, attr(TVP.XExo.filt, "X.tlag"), B)
        res.TVP.XExo[length(TVP.XExo.par.fin) + 5, 1] <- MSE.function(TVP.XExo.par.fin, 
            y, attr(TVP.XExo.filt, "X.tlag"), B)
        res.TVP.XExo[length(TVP.XExo.par.fin) + 6, 1] <- MASE.function(TVP.XExo.par.fin, 
            y, attr(TVP.XExo.filt, "X.tlag"), B)
        res.TVP.XExo[length(TVP.XExo.par.fin) + 7, 1] <- MSSE.function(TVP.XExo.par.fin, 
            y, attr(TVP.XExo.filt, "X.tlag"), B)
        colnames(res.TVP.XExo) <- c("Est.", "SE", "T test")
        rownames(res.TVP.XExo) <- c("mu0", "mu1", "sigma2", "pi11", 
            "pi22", "A1", "A2", "logLik", "AICc", "BIC", "MAE", 
            "MSE", "MASE", "MSSE")
        round(res.TVP.XExo, 3)
      
      
        TVP.XExo.smooth.series <- TVP.XExo.smooth[, 2]
        TVP.XExo.filt.p11 <- attr(TVP.XExo.filt, "p11")
        TVP.XExo.filt.p22 <- attr(TVP.XExo.filt, "p22")
      }
      
      
      

#Comparison

par.true

round(res.MS, 3)
round(res.TVP, 3)
round(res.GAS, 3)

if (type == "TVP-AR") round(res.TVP.XExo, 3)


ts.plot(cbind(pi11.true,MS.filt.p11,TVP.filt.p11,GAS.filt.p11)[(B+1):(B+size),],col=1:4,ylim=c(0,1.1),main="p11")
legend("top", legend=c("True", "Const.","TVP","GAS"),col=1:4, lty=1,horiz=T)

if (type == "TVP-AR") {

ts.plot(cbind(pi11.true,MS.filt.p11,TVP.filt.p11,GAS.filt.p11, TVP.XExo.filt.p11)[(B+1):(B+size),],col=1:5,ylim=c(0,1.1),main="p11")
legend("top", legend=c("True,", "Const.","TVP","GAS","TVP-XExo"),col=1:5, lty=1,horiz=T,cex=0.8)}





ts.plot(cbind(pi22.true,MS.filt.p22,TVP.filt.p22,GAS.filt.p22)[(B+1):(B+size),],col=1:4,ylim=c(0,1.1),main="p22")
legend("top", legend=c("True", "Const.","TVP","GAS"),col=1:4, lty=1,horiz=T)

if (type == "TVP-AR") {

ts.plot(cbind(pi22.true,MS.filt.p22,TVP.filt.p22,GAS.filt.p22, TVP.XExo.filt.p22)[(B+1):(B+size),],col=1:5,ylim=c(0,1.1),main="p22")
legend("top", legend=c("True", "Const.","TVP","GAS","TVP-XExo"),col=1:5, lty=1,horiz=T,cex=0.8)}
      




ts.plot(cbind(S.true-1 ,MS.smooth.series,TVP.smooth.series,GAS.smooth.series )[(B+1):(B+size),],col=1:4,ylim=c(0,1.1),main="Smoothing State 1")
legend("top", legend=c("True", "Const.","TVP","GAS"),col=1:4, lty=1,horiz=T)

if (type == "TVP-AR") {

ts.plot(cbind(S.true-1 ,MS.smooth.series,TVP.smooth.series,GAS.smooth.series, TVP.XExo.smooth.series)[(B+1):(B+size),],col=1:5,ylim=c(0,1.1),main="Smoothing State 1")
legend("top", legend=c("True", "Const.","TVP","GAS","TVP-XExo"),col=1:5, lty=1,horiz=T,cex=0.8)}


   
   
    
 


