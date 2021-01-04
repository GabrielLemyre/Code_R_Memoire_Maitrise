# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Kim Filtering
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
# KIM FILTER - ALGORITHME DE LISSAGE Implementation du filtre de Kim
#    Permet le calcul des probabilités lisées : phi_t = P[C_t | X_1:T]
#
# Either the filtered probs are given in the function call or the
# parameters mu and sigma required to calculate these probabilities
# --------------------------------------------------------
normal.HMM.KimFilter = function(DATA,
                                exp.var,
                                type,
                                distribution="Normal",
                                mu=NULL,
                                sigma=NULL,
                                matGamma,
                                p.ct.x1t.input=NULL,
                                initial.Distribution=NULL,
                                Transition.Type="Homogeneous",
                                nbStepsBack=0,
                                nu=NULL){
  diag.Gamma <- NULL
  
  # Testing if the parameters are given or if
  # the chain of filtered probabilities is given
  # Gamma must be specified either way
  if ((is.null(mu) || is.null(sigma)) && is.null(p.ct.x1t.input)){
    stop(cat("KIM FILTER : mu =",mu,'\nsigma =',sigma,'\np.ct.x1t =',p.ct.x1t.input,'\n'))
  }
  
  if (type %in% c("HHMM","HHMM.simplifie","DDMS")){
    Gamma.Built.obj <- Gamma.Build(prob.i=matGamma,
                                  type=type,
                                  Transition.Type=Transition.Type,
                                  nbStepsBack=nbStepsBack)
    
    Gamma.Built <- Gamma.Built.obj$matGamma
    
    if (type %in% c("DDMS")){
      mu <- mu %x% rep(1,nbStepsBack)
      sigma <- sigma %x% rep(1,nbStepsBack)
      initial.Distribution <- initial.Distribution %x% c(1,rep(0,nbStepsBack-1)) # Kronecker product
    }
  } else {
    Gamma.Built <- matGamma
  }
  
  # cat("type =",type,"-- Transition.Type =",Transition.Type,"-- dim(Gamma.Built) =",dim(Gamma.Built),"\n")
  # cat("type =",type,"-- Transition.Type =",Transition.Type,"-- length(Gamma.Built) =",length(Gamma.Built),"\n")
  
  n <- length(DATA)
  nbRegime <- length(mu)
  
  if (Transition.Type %in% c("LOGIT")){
    if (is.null(initial.Distribution)){
      stop(paste("KIM FILTER : The initial distribution must be provided explicitly if the Transitions are non-Homogeneous.\nPlease specify a vector of length ",
                 nbRegime," for the 'initial.Distribution' parameter."))
    }
  }
  
  # If the initial distribution P[C_1] isn't specified
  # the stationary distribution of Gamma is used
  if (is.null(initial.Distribution) & !type %in% c("UNIVARIEE","MELANGE") & !Transition.Type %in% c("GAS")){
    cat("solving for distr. init in Kim FILTER\n")
    initial.Distribution <- solve(t(diag(nbRegime)-Gamma.Built+1),rep(1,nbRegime))
  }
  
  
  if (type %in% c("MELANGE")){
    initial.Distribution <- Gamma.Built
  }
  
  if (type %in% c("UNIVARIEE")){
    initial.Distribution <- 1
  }
  
  # Calcul des probabilités filtrées par l'algorithme du filtre d'Hamilton
  # et retour de la diagonale de la matrice de transition si les transitions
  # ne sont pas homogènes
  Hamilton.Filter <- normal.HMM.HamiltonFilter(mu=mu,
                                               sigma=sigma,
                                               matGamma=Gamma.Built,
                                               initial.Distribution=initial.Distribution,
                                               DATA=DATA,
                                               exp.var = exp.var,
                                               distribution=distribution,
                                               Transition.Type=Transition.Type,
                                               nbStepsBack=nbStepsBack,
                                               nu=nu,
                                               type=type)
  
  
  if (!type %in% c("UNIVARIEE","MELANGE")){
    if (!Transition.Type %in% c("GAS")){
      
      # We keep the filtered probabilities from the Hamilton filter
      p.ct.x1t <- Hamilton.Filter$p.ct.x1t
      p.ct.x1tm1 <- Hamilton.Filter$p.ct.x1tm1
      
      # Initialization of the matrix to keep the diagonal of the transition matrix
      # step-by-step if the transitions aren't homogeneous
      diag.Gamma <- Hamilton.Filter$diag.Gamma
      
      if (Transition.Type!="Homogeneous"){
        if (sum(is.na(diag.Gamma))>0 & !Transition.Type %in% c("LOGIT")){
          err.message.builder(
            "KIM FILTER : NAs in diag.Gamma are at in Kim filter :\n",
            paste("1     ->",paste(round(diag.Gamma[,1],2),collapse=", "),"\n",
                  "(n-1) ->",paste(round(diag.Gamma[,(n-1)],2),collapse=", "),"\n",
                  "n     ->",paste(round(diag.Gamma[,n],2),collapse=", "),"\n")
          )
        }
      }
      
      # Alert message and information on localization of mistake if NAs are returned
      # by the Hamilton filter
      if (sum(is.na(p.ct.x1t))>0){
        err.message.builder(
          "KIM FILTER : NAs in p.ct.x1t are at in Kim filter :\n",
          paste("1     ->",paste(round(p.ct.x1t[,1],2),collapse=", "),"\n",
                "(n-1) ->",paste(round(p.ct.x1t[,(n-1)],2),collapse=", "),"\n",
                "n     ->",paste(round(p.ct.x1t[,n],2),collapse=", "),"\n")
        )
      }
      
      # Initialisation de la matrice de probabilités lissées
      p.ct.x1T <- matrix(0,nrow = nbRegime, ncol = n)
      
      # Probabilite de la derniere observation de la variable latente
      p.ct.x1T[,n] <- p.ct.x1t[,n]
      if (sum(is.na(p.ct.x1t[,n]))>0){
        stop("Mauvaise initialisation de p.ct.x1t[,n] =",p.ct.x1t[,n],"\n")
      }
    
      if (Transition.Type %in% c("LOGIT","LOGIT.w.filter")){
        if (Transition.Type=="LOGIT"){
          exp.var.Gamma <- if (is.null(exp.var)){
            DATA[n]
          } else {
            exp.var[n]
          }
        } else if (Transition.Type=="LOGIT.w.filter"){
          exp.var.Gamma <-if (is.null(exp.var)){
            c(DATA[n],p.ct.x1t[,n])
          } else {
            c(exp.var[n],p.ct.x1t[,n])
          }
        }
        Gamma.Built.obj <- Gamma.Build(prob.i = matGamma,
                                       Transition.Type=Transition.Type,
                                       nbStepsBack=nbStepsBack,
                                       exp.var=exp.var.Gamma,
                                       type=type)
        
        diag.Gamma[,n] <- Gamma.Built.obj$diag.Gamma
        Gamma.Built <- Gamma.Built.obj$matGamma
      }
      
      # Algorithme de lissage (Filtre de Kim)
      for (t in (n-1):1){
        # phi.t.ij <- matrix(NA,nbRegime,nbRegime)
        if (Transition.Type %in% c("LOGIT","LOGIT.w.filter")){
          if (Transition.Type=="LOGIT"){
            exp.var.Gamma <- if (is.null(exp.var)){
              DATA[t]
            } else {
              exp.var[t]
            }
          } else if (Transition.Type=="LOGIT.w.filter"){
            exp.var.Gamma <-if (is.null(exp.var)){
              c(DATA[t],p.ct.x1t[,t])
            } else {
              c(exp.var[t],p.ct.x1t[,t])
            }
          }
          Gamma.Built.obj <- Gamma.Build(prob.i = matGamma,
                                         Transition.Type=Transition.Type,
                                         nbStepsBack=nbStepsBack,
                                         exp.var=exp.var.Gamma,
                                         type=type)
          
          diag.Gamma[,t] <- Gamma.Built.obj$diag.Gamma
          Gamma.Built <- Gamma.Built.obj$matGamma
        }
        
        # cat("In kim step at time", t,": dim(Gamma.Built) =",dim(Gamma.Built),"\n")
        
        # Version VECTORIELLE
        # p.ct.x1T[,t] <- as.numeric(p.ct.x1t[,t]*((p.ct.x1T[,t+1]/p.ct.x1t[,t+1])%*%t(Gamma.Built)))
        p.ct.x1T[, t] <- (Gamma.Built %*% (p.ct.x1T[,t+1]/p.ct.x1tm1[,t+1])) * p.ct.x1t[,t]
        
      }
      
    } else { # Si modèle GAS
      # Extraction info du filtre d'Hamilton
      p.ct.x1t <- Hamilton.Filter$p.ct.x1t
      p.ct.x1tm1 <- Hamilton.Filter$p.ct.x1tm1
      diag.Gamma <- Hamilton.Filter$diag.Gamma
      
      # Obtention des composantes de Gamma
      p11 <- diag.Gamma[1,]
      p22 <- diag.Gamma[2,]
      
      p12 <- 1 - p11
      p21 <- 1 - p22
      
      p.ct.x1T <- matrix(0, ncol = n, nrow = 2)
      p.ct.x1T[, n] <- p.ct.x1t[, n]
      
      for (t in ((n - 1):1)) {
        Gamma.Built <- matrix(c(p11[t + 1], p12[t + 1], 
                                p21[t + 1], p22[t + 1]), 
                              ncol = 2, 
                              byrow = T)
        # p.ct.x1T[, t] <- p.ct.x1t[,t]*((p.ct.x1T[,t+1]/(p.ct.x1t[,t+1]))%*%t(Gamma.Built))
        
        p.ct.x1T[, t] <- (Gamma.Built %*% (p.ct.x1T[,t+1]/p.ct.x1tm1[,t+1])) * p.ct.x1t[,t]
      }
      
    }
  } else { # Si UNIVARIEE OU MELANGE
    
    # We keep the filtered probabilities from the Hamilton filter
    p.ct.x1t <- Hamilton.Filter$p.ct.x1t
    
    # Initialization of the matrix to keep the diagonal of the transition matrix
    # step-by-step if the transitions aren't homogeneous
    diag.Gamma <- Hamilton.Filter$diag.Gamma
    
    # Pas de changement de régimes
    p.ct.x1T <- Hamilton.Filter$p.ct.x1t
    # print(Hamilton.Filter$p.ct.x1t)
  }
  
  # Alert message and information on localization of mistake if NAs are returned
  # by the Hamilton filter
  if (sum(is.na(p.ct.x1T))>0){
    err.message.builder(
      "KIM FILTER : NAs in p.ct.x1T are at in Kim filter :\n",
      paste("1     ->",paste(round(p.ct.x1T[,1],2),collapse=", "),"\n",
            "(n-1) ->",paste(round(p.ct.x1T[,(n-1)],2),collapse=", "),"\n",
            "n     ->",paste(round(p.ct.x1T[,n],2),collapse=", "),"\n")
    )
  }
  
  return(list(llk=Hamilton.Filter$llk,
              p.ct.x1T=p.ct.x1T,
              p.ct.x1t=Hamilton.Filter$p.ct.x1t,
              p.ct.x1tm1=Hamilton.Filter$p.ct.x1tm1,
              diag.Gamma = diag.Gamma,
              u.t=Hamilton.Filter$u.t,
              esp.xt=Hamilton.Filter$esp.xt,
              pred.err.xt=Hamilton.Filter$pred.err.xt,
              standard.pred.err.xt=Hamilton.Filter$standard.pred.err.xt,
              diag.Gamma=Hamilton.Filter$diag.Gamma,
              initial.Distribution = Hamilton.Filter$initial.Distribution)
  )
  
}


