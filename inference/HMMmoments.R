# ----------------------------------------------------------------------------------------------------
# HIDDEN MARKOV MODELS
# MOMENTS (MEAN, VARIANCE, AUTOCOVARIANCE AND AUTOCORRELATION)
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# First version : 28 février 2020
# Last version  : 28 février 2020
# ----------------------------------------------------------------------------------------------------

VarianceHMM <- function(Type.Transition="Homogeneous", dist.initiale, ...){
  opt.args <- list(...)
  
  # Unpacking optional arguments --------------------------------------------------------------------------
  if (length(opt.args)>0){
    for(i in 1:length(opt.args)) {
      assign(x = names(opt.args)[i], value = opt.args[[i]])
    }
  }
  
  
  # Definition  --------------------------------------------------------------------------
  K <- dim(Gamma)[1] # Dimension de l'espace d'état
  
  upsilon <- 1:K
  Upsilon <- diag(upsilon)
  
  vecUn <- rep(1,K)
  
  if (Transition.Type=="Homogeneous"){
    distribution.Distribution <- vecUn%*%solve(diag(K)-Gamma+(vecUn%*%t(vecUn))) 
    
    if(is.null(initial.Distribution)){
      
      # Utilisation de la distribution stationnaire comme distribution initiale
      initial.Distribution <- distribution.Distribution
      
      if (nbStepsEsp>1){
        cat("Avertissement : \n",
            "  Puisque la distribution stationnaire est utilisée comme \n",
            "  distribution initiale et puisque celle-ci multipliée par \n",
            "  la matrice de transition donne elle même, l'espérance de \n",
            "  chaque pas est la même et ne sera donc calculée qu'une fois.\n",
            "  Ainsi la variable 'nbStepsEsp=",nbStepsEsp,"' est remplacée par 1.\n",sep="")
        
        nbStepsEsp <- 1
      }
    }
  }
  
}