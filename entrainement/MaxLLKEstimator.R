# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Maximum Likelihood Estimator
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# Last version : january 9th, 2020
# Last version : january 9th, 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# OPTIMISATION Fonction afin d'optimiser numériquement la log-vraisemblance
#   Optimisation à l'aide du filtre d'Hamilton (Probabilités filtrées
#   P[C_t | X_1:t])
# --------------------------------------------------------
normal.HMM.mle <- function(data.var,
                           exp.var,
                           mu, 
                           sigma, 
                           matGamma, 
                           type=NULL, 
                           distribution,
                           Transition.Type,
                           nbStepsBack, nu,
                           initial.Distribution,
                           print.result.optim=TRUE){
  nbRegime <- length(mu)
  parvect0       <- normal.HMM.N2W(mu, sigma, matGamma, 
                                   type=type,
                                   Transition.Type=Transition.Type)
  
  # cat("length(data.var) =",length(data.var),"\n")
  # cat("length(exp.var) =",length(exp.var),"\n")
  
  #Calcul du EVM
  start_time <- Sys.time()
  mod            <- optim(par=parvect0,
                          fn=normal.HMM.mllk,
                          method="BFGS",
                          DATA=data.var,
                          exp.var=exp.var,
                          type=type,
                          distribution=distribution,
                          nu=nu,
                          Transition.Type=Transition.Type,
                          nbStepsBack=nbStepsBack,
                          initial.Distribution=initial.Distribution,
                          control=list(maxit=1e+06,
                                       reltol=1e-8), 
                          hessian = TRUE)
  
  se.vec <- NULL
  tryCatch({
    cat("Diagonale de la matrice hessienne :\n")
    print(mod$hessian)
    vc.mat <- solve(mod$hessian) # var-cov matrix
    se.vec <- sqrt(diag(vc.mat))        # standard errors
  },
  error = function(err){
    err.message.builder("La matrice hessienne ne semble pas inversible :\n",
                        err)
  },
  warning = function(warn.mess){
    warn.message.builder("La matrice hessienne ne semble pas inversible :\n",
                         warn.mess)
  })
  
  # fisher_info <- solve(fit$hessian)
  # prop_sigma <- sqrt(diag(fisher_info))
  # prop_sigma <- diag(prop_sigma)
  # upper <- mod$par + 1.96*prop_sigma
  # lower <- mod$par - 1.96*prop_sigma
  # interval <- data.frame(value=mod$par, upper=upper, lower=lower)
  
  # mod            <- nlm(normal.HMM.mllk,
  #                       parvect0,
  #                       DATA=data.var,
  #                       type=type,
  #                       distribution=distribution,
  #                       nu=nu,
  #                       Transition.Type=Transition.Type,
  #                       nbStepsBack=nbStepsBack,
  #                       initial.Distribution=initial.Distribution,
  #                       iterlim=1000,
  #                       gradtol = 1e-10)
  # mod            <- nlminb(start=parvect0,
  #                         objective=normal.HMM.mllk,
  #                         DATA=data.var,
  #                         type=type,
  #                         distribution=distribution,
  #                         nu=nu,
  #                         Transition.Type=Transition.Type,
  #                         nbStepsBack=nbStepsBack,
  #                         initial.Distribution=initial.Distribution,
  #                         control=list(eval.max=500,
  #                                      iter.max=300,
  #                                      rel.tol=1e-10))
  
  
  # natural.par.optim <- normal.HMM.W2N(mod$estimate,
  #                                     type,
  #                                     Transition.Type=Transition.Type) # [NLM]
  natural.par.optim <- normal.HMM.W2N(mod$par,
                                      type,
                                      Transition.Type=Transition.Type) # [NLMINB or OPTIM]
  
  if (distribution=='Student'){
    optim.sigma=natural.par.optim$sigma*sqrt(nu/(nu-2))
  } else {
    optim.sigma=natural.par.optim$sigma
  }
  
  
  if (print.result.optim){
    print(mod)
    
    end_time <- Sys.time()
    hline()
    cat("Time for optimization :",end_time - start_time,"\n")
  }
  
  #Log-vraisemblance
  # mllk           <- -mod$minimum # Puisqu'on minimise -llk, on veut tout de même retourner la valeur à être maximisée [NLM]
  mllk           <- -mod$value # Puisqu'on minimise -llk, on veut tout de même retourner la valeur à être maximisée [OPTIM]
  # mllk           <- -mod$objective # Puisqu'on minimise -llk, on veut tout de même retourner la valeur à être maximisée [NLMINB]
  np             <- length(parvect0)
  AIC.v            <- 2*(-mllk+np)
  n              <- sum(!is.na(data.var))
  AICc           <- AIC.v + (2*np^2 + 2*np)/(n-np-1)
  BIC.v            <- np*log(n)-2*mllk
  # codeConv           <- mod$code # [NLM]
  codeConv           <- mod$convergence # [OPTIM, NLMINB]
  messageConv        <- mod$message # [OPTIM, NLMINB]
  countsConv        <- mod$counts # [OPTIM, NLMINB]
  MAFE=0
  MSFE=0
  MASFE=0
  MSSFE=0
  
  return(list(
    nb.obs=n,
    mu=natural.par.optim$mu, 
    sigma=optim.sigma, 
    matGamma=natural.par.optim$matGamma, 
    nb.param = length(parvect0),
    mllk=mllk,
    AIC.v=AIC.v, 
    AICc=AICc, 
    BIC.v=BIC.v,
    MAFE=MAFE,
    MSFE=MSFE,
    MASFE=MASFE,
    MSSFE=MSSFE,
    codeConv=codeConv,
    messageConv=messageConv,
    countsConv=countsConv,
    mod=mod,
    se.vec=se.vec
  ))
}