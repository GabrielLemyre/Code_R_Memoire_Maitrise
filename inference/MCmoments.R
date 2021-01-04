# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Moments of the Markov Chain
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# Last version : january 10th, 2020
# Last version : january 13th, 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# MOMENTS DE LA CHAÎNE DE MARKOV
# --------------------------------------------------------
# Fonction pour calculer les moments des v.a. de la chaîne de Markov
#     Calcul de l'espérance, de l'autocovariance et autocorrélation 
#     de chaîne homogène/stationnaire ou non.
# --------------------------------------------------------
MomentsMC = function(Gamma.matrix,                         # Matrice de transition de la chaîne de Markov
                     initial.Distribution=NULL,     # Distribution initiale de la MC. Si absente, utilisation de la dist asymptotique
                     Transition.Type="Homogeneous", # Type de transition
                     Data=NULL,                     # Si transition non-homogène, nécessaire pour construire matrice de transition
                     nbStepsEsp=100,                # Nombre de pas en avant pour calcul de l'espérance dans la MC
                     nbStepsACF=20){                # Nombre de pas pour le calcul de Cov[C_t, C_(t+k)]
  
  K <- dim(Gamma.matrix)[1] # Dimension de l'espace d'état
  
  upsilon <- 1:K
  Upsilon <- diag(upsilon)
  
  vecUn <- rep(1,K)
  
  if (Transition.Type=="Homogeneous"){
    distribution.stationnaire <- vecUn%*%solve(diag(K)-Gamma.matrix+(vecUn%*%t(vecUn))) 
    
    if(is.null(initial.Distribution)){
      
      # Utilisation de la distribution stationnaire comme distribution initiale
      initial.Distribution <- distribution.stationnaire
      
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
  
  Esperance.Ct <- rep(NA,nbStepsEsp) # Initialisation du vecteur des espérances
  ACF.Ct <- rep(NA,nbStepsACF) # Initialisation du vecteur des espérances
  Gamma.matrix.m.step <- diag(K) # Matrice identité pour premier pas
  
  for (i in 1:max(nbStepsEsp,nbStepsACF)){
    # E[Ct] = u(1) %*% Gamma.matrix^(t-1) %*% c(1, ..., K)
    Gamma.matrix.m.step <- Gamma.matrix.m.step%*%Gamma.matrix
    if (i<=nbStepsEsp){
      Esperance.Ct[i] <- initial.Distribution%*%Gamma.matrix.m.step%*%upsilon
    }
    if (i<=nbStepsACF){
      ACF.Ct[i] <- initial.Distribution%*%Gamma.matrix.m.step%*%Upsilon%*%Gamma.matrix.m.step%*%upsilon-(initial.Distribution%*%Gamma.matrix.m.step%*%upsilon)^2
    }
    
  }
  
  # plot(Esperance.Ct, ylim=c(0.999,K+0.001), type="l")
  # plot(ACF.Ct, type="h")
  
  Esperance.Ct.df <- data.frame(c(1:length(Esperance.Ct)), Esperance.Ct)
  names(Esperance.Ct.df) <- c("t","esp")
  
  ACF.Ct.df <- data.frame(c(1:length(ACF.Ct)), Esperance.Ct)
  names(Esperance.Ct.df) <- c("t","esp")
  
  ACF.Ct.df <- with(ACF.Ct, data.frame(lag, acf))
  
  Esperance.plot <- ggplot(data = Esperance.Ct.df, mapping = aes(x = t, y = esp)) +
    geom_segment(mapping = aes(xend = t, yend = 0)) +
    labs(title="Différence entre Corrélation valeur absolue et au carré",
         x ="Lags", 
         y = "Différence") + 
    theme(
      plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
      axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
      axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
    )
  
  ACov.plot <- ggplot(data = ACF.Ct.df, mapping = aes(x = lag, y = acf)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) +
    labs(title="Différence entre Corrélation valeur absolue et au carré",
         x ="Lags", 
         y = "Différence") + 
    theme(
      plot.title = element_text(color="Black", size=title.size.var, face="bold.italic"),
      axis.title.x = element_text(color="black", size=axis.lab.size.var, face="bold"),
      axis.title.y = element_text(color="black", size=axis.lab.size.var, face="bold")
    )
  
  p <- plot_grid(Esperance.plot,ACov.plot, labels = "AUTO", ncol = 1)
  # ggsave(paste(GraphPath,"/1 - DataStats/Dataset_ACF_",index,".png",sep=""), p)
  
  MomentsMC <- list(Esperance.Ct=Esperance.Ct,
                    ACF.Ct=ACF.Ct,
                    distribution.stationnaire=distribution.stationnaire)
}
# test