# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Unconstraining parameters
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
# PARAMÈTRES CONTRAINTS -> PARAMÈTRES NON-CONTRAINTS
# --------------------------------------------------------
# Function to convert from natural parameters to their working counterparts
#     The idea is that by changing the constrained parameters to
#     unconstrained versions, the optimization can be done without
#     constraints
# --------------------------------------------------------
normal.HMM.N2W = function(mu, sigma, matGamma=NULL,
                          type, Transition.Type="Homogeneous"){
  
  if (length(mu)!=length(sigma)){
    stop("Les dimensions des vecteurs mu et sigma ne concordent pas.")
  }
  
  if(!type %in% c("UNIVARIEE","HHMM","HHMM.simplifie", "DDMS")){
    if (!Transition.Type %in% c("LOGIT","GAS")){
      dim.gamma <- if (!is.null(dim(matGamma)[1])){
        dim(matGamma)[1]
      } else {
          length(matGamma)
        }
      if (length(mu)!=dim.gamma){
        stop("Les dimensions de la matrice matGamma et des vecteurs mu et sigma ne concordent pas.")
      }
    }
  }
  
  tsigma <- log(sigma)
  tmu <- mu # Aucune transformation nécessaire
  tGamma <- NULL # Construit plus loin
  
  if (type=='HMM'){
    if ((Transition.Type=="Homogeneous")||is.null(Transition.Type)){
      # --------------
      # CLASSIC HMM
      # --------------
      nbRegime <- length(mu)
      trans.Gamma <- log(matGamma/matGamma[,nbRegime])
      # Retrait de la dernière colonne de notre matrice de transition
      tGamma      <- as.vector(t(trans.Gamma[,-nbRegime]))
    } else if (Transition.Type=="LOGIT" || Transition.Type=="LOGIT.w.filter"){
      # --------------
      # NON-HOMOGENOUS
      # --------------
      nbRegime <- length(mu)
      tGamma <- matGamma # Because Beta_0 and Beta_1 are unconstrained vectors of parameters
    } else if (Transition.Type=="GAS"){
      # --------------
      # GAS -- NON-HOMOGENOUS
      # --------------
      nbRegime <- length(mu)
      tGamma <- log(matGamma) - log(1 - matGamma) # 
    }
    
  } else if (type=="DDMS"){
    # --------------
    # NON-HOMOGENOUS
    # --------------
    nbRegime <- length(mu)
    tGamma <- matGamma # Because Beta_0 and Beta_1 are unconstrained vectors of parameters
  } else if (type=='HHMM'){
    # ----------------
    # HIERARCHICAL HMM
    # ----------------
    Gamma.temp <- matGamma[1:4,]
    trans.Gamma <- log(Gamma.temp/Gamma.temp[,3])
    tGamma.temp <- trans.Gamma[,1:2]
    
    for (i in 1:4){
      tGamma=cbind(tGamma, t(tGamma.temp[i,]))
    } 
    
    Gamma.temp <- matGamma[5:8,1:2]
    trans.Gamma <- log(Gamma.temp/Gamma.temp[,2])
    tGamma.temp <- matrix(trans.Gamma[,1],ncol=1)
    
    for (i in 1:4){
      tGamma=cbind(tGamma, tGamma.temp[i,])
    }
  } else if (type=='HHMM.simplifie'){
    nbRegime <- length(mu)
    trans.Gamma <- log(matGamma/matGamma[,2])
    # Retrait de la dernière colonne de notre matrice de transition
    tGamma      <- as.vector(t(trans.Gamma[,-2]))
  } else if (type=='MELANGE'){
    nbRegime <- length(mu)
    trans.Gamma <- log(matGamma/matGamma[nbRegime])
    # Retrait de la dernière colonne de notre matrice de transition
    tGamma      <- trans.Gamma[-nbRegime]
  }
  
  parvect = c(tmu,tsigma)
  if(!is.null(tGamma)){
    parvect <- c(parvect,tGamma)
  }
  
  return(parvect)
}
