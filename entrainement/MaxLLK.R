# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Maximum Likelihood computation
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
# FONCTION OBJECTIF afin de calculer la vraisemblance et la log-vraisemblance
#     Le calcul est effectué à l'aide du filtre d'hamilton (Probabilités
#     filtrées P[C_t | X_1:t])
# --------------------------------------------------------
normal.HMM.mllk = function(parvect, 
                           DATA, exp.var,
                           type, distribution="Normal",Transition.Type="Homogeneous",
                           nbStepsBack=0, nu=NULL,
                           initial.Distribution){
  
  n <- length(DATA)
  
  natural.par <- normal.HMM.W2N(parvect=parvect,
                                type=type,
                                Transition.Type=Transition.Type)
  
  # Assignation des paramètres naturels
  mu <- natural.par$mu
  sigma <- natural.par$sigma
  matGamma <- natural.par$matGamma
  
  if (!type %in% c("UNIVARIEE")){
    
    if (type %in% c("HHMM","HHMM.simplifie")){
      Gamma.vec <- Gamma.Build(prob.i=matGamma,
                               type=type,
                               Transition.Type=Transition.Type,
                               nbStepsBack=nbStepsBack)
      
      matGamma <- Gamma.vec$matGamma
    }
    
    if (type=="DDMS"){
      Gamma.vec <- Gamma.Build(prob.i=matGamma,
                               type=type,
                               Transition.Type=Transition.Type,
                               nbStepsBack=nbStepsBack)
      
      matGamma <- Gamma.vec$matGamma
      
      mu <- mu %x% rep(1,nbStepsBack)
      sigma <- sigma %x% rep(1,nbStepsBack)
      initial.Distribution <- initial.Distribution %x% c(1,rep(0,nbStepsBack-1)) # Kronecker product
    }
    
    optim.results <- normal.HMM.HamiltonFilter(mu=mu,
                                               sigma=sigma,
                                               matGamma=matGamma,
                                               DATA=DATA,
                                               exp.var=exp.var,
                                               distribution=distribution,
                                               nu=nu,
                                               Transition.Type=Transition.Type,
                                               type=type,
                                               nbStepsBack=nbStepsBack,
                                               initial.Distribution=initial.Distribution)
    # print(optim.results$llk)
    llk <- -optim.results$llk
  } else {
    llk <- -sum(dnorm(DATA, 
                      mean = mu, 
                      sd = sigma, 
                      log = TRUE))
  }
  return(llk)
}