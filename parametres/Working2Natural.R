# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Constraining parameters
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# Last version : january 8th, 2020
# Last version : january 9th, 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# PARAMÈTRES NON-CONTRAINTS -> PARAMÈTRES CONTRAINTS
# --------------------------------------------------------
# Function to convert from Working parameters to their natural equivalents
# --------------------------------------------------------
normal.HMM.W2N = function(parvect, type, Transition.Type){
  # On initialise la matrice Gamma
  matGamma <- NULL
  
  if (type=="HMM"){
    if ((Transition.Type=="Homogeneous")||is.null(Transition.Type)){
      # --------------
      # CLASSIC HMM
      # --------------
      nparvect=length(parvect)
      nbRegime=(-1+sqrt(1+4*nparvect))/2
      
      nat.param  <- exp(parvect)
      mu         <- log(nat.param[1:nbRegime])
      sigma      <- nat.param[(nbRegime+1):(2*nbRegime)]
      
      matGamma   <- t(matrix(nat.param[(2*nbRegime+1):(2*nbRegime+nbRegime*(nbRegime-1))],
                          nrow = nbRegime-1, ncol = nbRegime))
      matGamma   <- cbind(matGamma, rep(1, nbRegime))
      matGamma   <- matGamma/apply(matGamma,1,sum)
      
    } else if (Transition.Type=="LOGIT" || Transition.Type=="LOGIT.w.filter"){
      # --------------
      # NON-HOMOGENOUS
      # --------------
      nparvect=length(parvect)
      nbRegime=4
      
      nat.param  <- exp(parvect)
      mu         <- log(nat.param[1:nbRegime])
      sigma      <- nat.param[(nbRegime+1):(2*nbRegime)]
      
      matGamma   <- log(nat.param[(2*nbRegime+1):length(nat.param)])
    } else if (Transition.Type=="GAS"){
      # --------------
      # GAS -- NON-HOMOGENOUS
      # --------------
      nparvect=length(parvect)
      nbRegime=2
      
      nat.param  <- parvect
      mu         <- nat.param[1:nbRegime]
      sigma      <- exp(nat.param[(nbRegime+1):(2*nbRegime)])
      
      gamma.terms <- nat.param[(2*nbRegime+1):length(nat.param)]
      matGamma   <- exp(gamma.terms)/(1 + exp(gamma.terms))
    }
    
  } else if (type=="DDMS"){
    # --------------
    # DURATION-DEPENDANT-MARKOV-SWITCHING
    # --------------
    nparvect=length(parvect)
    nbRegime=4
    
    nat.param  <- exp(parvect)
    mu         <- log(nat.param[1:nbRegime])
    sigma      <- nat.param[(nbRegime+1):(2*nbRegime)]
    
    matGamma   <- log(nat.param[(2*nbRegime+1):length(nat.param)])
  } else if (type=="FHMM"){ # TODO
    # --------------
    # FACTORIAL HMM
    # --------------
    nbRegime <- 4
    mu <- parvect[1:4]
    
    nat.param  <- exp(parvect)
    sigma <- nat.param[5:8]
    
    matGamma <- matrix(nat.param[9:12],ncol=2,byrow=T)
    matGamma <- cbind(matGamma, rep(1, nbRegime))
    Gamma.temp.1 <- matGamma/apply(matGamma,1,sum)
    
    
  } else if (type=="HHMM"){ 
    # --------------
    # HIERARCHICAL HMM
    # --------------
    nbRegime <- 4;
    mu <- parvect[1:4];
    
    nat.param  <- exp(parvect)
    sigma <- nat.param[5:8]
    
    matGamma <- matrix(nat.param[9:16],ncol=2,byrow=T)
    matGamma <- cbind(matGamma, rep(1, nbRegime))
    Gamma.temp.1 <- matGamma/apply(matGamma,1,sum)
    
    matGamma <- matrix(nat.param[17:20],ncol=1,byrow=T)
    matGamma <- cbind(matGamma, rep(1, nbRegime))
    Gamma.temp.2 <- matGamma/apply(matGamma,1,sum)
    
    # Padding the second part to fit dimensions of the first one and
    # then stacking them
    matGamma = rbind(Gamma.temp.1, cbind(Gamma.temp.2,rep(0, nbRegime)));
  } else if (type =="HHMM.simplifie"){
    nparvect=length(parvect)
    nbRegime=4
    
    nat.param  <- exp(parvect)
    mu         <- log(nat.param[1:nbRegime])
    sigma      <- nat.param[(nbRegime+1):(2*nbRegime)]
    
    matGamma   <- matrix(nat.param[(2*nbRegime+1):nparvect],
                           ncol = 1)
    
    matGamma   <- cbind(matGamma, rep(1, 6))
    matGamma   <- matGamma/apply(matGamma,1,sum)
  } else {
    nat.param  <- exp(parvect)
    if (type=="MELANGE"){
      nbRegime=(length(nat.param)+1)/3
      matGamma <- c(nat.param[(2*nbRegime+1):length(nat.param)],1)
      matGamma   <- matGamma/sum(matGamma)
    } else {
      nbRegime=1
    }
    mu         <- log(nat.param[1:nbRegime])
    sigma      <- nat.param[(nbRegime+1):(2*nbRegime)]
  }
  
  return(list(mu=mu, sigma=sigma, matGamma=matGamma))
}