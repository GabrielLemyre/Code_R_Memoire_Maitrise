# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Training for the parameters
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# Last version : 15 avril 2019
# Last version : 4 mars 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# FULL HMM EXPERIENCE : TRAINING FOR THE PARAMETERS
# --------------------------------------------------------
# 
# --------------------------------------------------------
train.HMM = function(DATA,
                     data.name,
                     exp.var.name=NULL,
                     mu0=NULL,
                     sigma0=NULL,
                     gamma0=NULL,
                     distribution="Normal",
                     initial.Distribution=NULL,
                     nom.modele.row="",
                     type,
                     Transition.Type="Homogeneous",
                     nbStepsBack=1,                 # Seulement si non-homogène
                     nu=4,                          # Seulement si distribution student
                     nbRegime=NULL,                 # Seulement si auto-assign=TRUE
                     auto.assign=FALSE,
                     print.result.optim=TRUE){
  
  
  Observation.Driven.models <- c("LOGIT")
  hequal()
  cat("DÉBUT : ENTRAÎNEMENT pour le \n     [ ",
      nom.modele.row,
      " ]\n")
  hequal()
  
  if (type %in% c("HHMM","HHMM.simplifie")){
    nbRegime <- 4
  } else if (is.null(nbRegime)){
    nbRegime <- length(mu0)
  }
  
  # Matrice de coefficients pour LOGIT
  A <- matrix(nrow = nbRegime, ncol = nbRegime)
  B <- matrix(nrow = nbRegime, ncol = nbRegime)
  
  # fail-safe au cas où erreur lors de l'utilisation des modèles
  #   non-homogènes
  if (Transition.Type!="Homogeneous"){
    nbRegime <- length(mu0)
  }
  # --------------------------------------------
  # TRAINING FOR OPTIMAL PARAMETERS (HAMILTON)
  # --------------------------------------------
  
  start_time <- as.numeric(Sys.time()) # Start timer
  # start_time <- get.Sys.Time() # Start timer
  
  # Selection des données pour l'entrainement
  data.var <- as.numeric(as.character(unlist(DATA[data.name])))
  exp.var <- if(Transition.Type %in% Observation.Driven.models){
    if (!is.null(exp.var.name)){
      as.numeric(as.character(unlist(DATA[exp.var.name])) )
    } else {
      as.numeric(as.character(unlist(DATA[data.name])) )
    }
  }else{
    NULL 
  }
  
  if (!length(data.var)>0){stop("FATAL ----\n\n              Before train length(data.var) = ",length(data.var),"\ndata.name = ",data.name,"\n")}
  if (!length(exp.var)>0 & Transition.Type %in% Observation.Driven.models){stop("Des observations devraient se rendre à train length(exp.var) = ",length(exp.var))}
  if (length(exp.var)>0 & !Transition.Type %in% Observation.Driven.models){stop("Il ne devrait pas y avoir d'observations qui se rendent à train length(exp.var) = ",length(exp.var))}
  # stop("Before train nbRegime =",nbRegime,"\n")
  
  HMM.Train <- normal.HMM.mle(data.var = data.var,
                              exp.var = exp.var,
                              mu=mu0,
                              sigma=sigma0,
                              matGamma=gamma0,
                              type=type,
                              distribution=distribution,
                              nu=nu,
                              Transition.Type=Transition.Type,
                              nbStepsBack=nbStepsBack,
                              initial.Distribution=initial.Distribution,
                              print.result.optim=print.result.optim)
  
  # print("Post train")
  end_time <- as.numeric(Sys.time()) # End timer
  # end_time <- get.Sys.Time() # End timer
  
  timeSpent.String <- get.Time.Diff(start_time, end_time)
  
  # If the model is the basic HMM, the parameters are ordered in decreasing
  # order of the volatility values
  # if (type=="HMM" & !Transition.Type %in% c("LOGIT","GAS")){
  #   HMM.Train <- order.sigma(HMM.Train,Transition.Type=Transition.Type)
  # }
  
  print(HMM.Train[!names(HMM.Train) %in% c("MAFE","MSFE","MASFE","MSSFE")])
  
  # --------------------------------------------
  # KIM FILTER -> SMOOTHED PROBABILITIES P[C_t|x_1:T]
  # --------------------------------------------
  # print("Post Train")
  # if (!type %in% c("UNIVARIEE","MELANGE")){
  # cat("Longueur de data.var avant kim length(data.var) =",length(data.var),"\n")
  smooth.Prob <- normal.HMM.KimFilter(DATA = data.var,
                                      exp.var = exp.var,
                                      type=type,
                                      distribution=distribution,
                                      mu=HMM.Train$mu,
                                      sigma=HMM.Train$sigma,
                                      matGamma=HMM.Train$matGamma,
                                      nu=nu,
                                      Transition.Type=Transition.Type,
                                      nbStepsBack=nbStepsBack,
                                      initial.Distribution=initial.Distribution)
  
  
  
  pred.err.xt=smooth.Prob$pred.err.xt
  standard.pred.err.xt=smooth.Prob$standard.pred.err.xt
  
  # ——————————————
  # IS - MAFE
  HMM.Train$MAFE <- 
    mean(abs(pred.err.xt))
  # ——————————————
  
  
  # ——————————————
  # IS - MSFE
  HMM.Train$MSFE <- 
    mean((pred.err.xt)^2)
  # ——————————————
  
  
  # ——————————————
  # IS - MASFE
  HMM.Train$MASFE <- 
    mean(abs(standard.pred.err.xt))
  # ——————————————
  
  
  # ——————————————
  # IS - MSSFE
  HMM.Train$MSSFE <- 
    mean((standard.pred.err.xt)^2)
  # ——————————————
  
  hline(20)
  cat("MAFE  : ",round(HMM.Train$MAFE,3),"\nMSFE  : ",round(HMM.Train$MSFE,3),"\n",
      "MASFE : ",round(HMM.Train$MASFE,3),"\nMSSFE : ",round(HMM.Train$MSSFE,3),"\n",sep="")
  
  
  if (type %in% c("HHMM","HHMM.simplifie","HMM","DDMS")){
    hline(60)
    cat("Matrice de transition résultante :\n")
    hline(60)
    
    if (type %in% c("HHMM","HHMM.simplifie")){
      matGamma <- Gamma.Build(HMM.Train$matGamma, 
                              type = type,
                              Transition.Type=Transition.Type)
      print(list(
        `Matrice de transition` = round(matGamma$matGamma,3),
        `Diagonale principale` = round(matGamma$diag.Gamma,3)
      ))
      
    } else if (type == "DDMS"){
      matGamma <- Gamma.Build(HMM.Train$matGamma, 
                              type = type,
                              Transition.Type=Transition.Type,
                              nbStepsBack = nbStepsBack)
      print(list(
        `Coefficients beta0` = round(matGamma$beta0,3),
        `Coefficients beta1` = round(matGamma$beta1,3)
      ))
      
    } else if (Transition.Type == "LOGIT"){
      matGamma <- Gamma.Build(HMM.Train$matGamma, 
                              type = type,
                              Transition.Type="LOGIT.final")
      print(list(
        `Coefficients beta0` = round(matGamma$beta0,3),
        `Coefficients beta1` = round(matGamma$beta1,3)
      ))
      
    } else if (Transition.Type == "GAS"){
      print(
        list(
          `Distribution initiale` = round(smooth.Prob$initial.Distribution,3),
          A = round(diag(HMM.Train$matGamma[3:4]),3),
          B = round(diag(HMM.Train$matGamma[5:6]),3)
        )
      )
    } else if (type == "HMM"){
      print(list(
        `Matrice de transition` = round(HMM.Train$matGamma,3),
        `Diagonale principale` = round(diag(HMM.Train$matGamma),3)
      ))
    }
  } 
  
  hline()
  cat("Temps d'entraînement :",timeSpent.String,"\n",sep="")
  hline()
  
  hequal()
  cat("     [ ",
      nom.modele.row,
      " ]\n")
  cat("FIN : ENTRAÎNEMENT\n")
  hequal()
  
  return(list(nom.modele.row=nom.modele.row,
              HMM.Train=HMM.Train,
              smooth.Prob=smooth.Prob$p.ct.x1T,
              filtered.Prob=smooth.Prob$p.ct.x1t,
              diag.Gamma=smooth.Prob$diag.Gamma,
              pred.err.xt=pred.err.xt,
              standard.pred.err.xt=standard.pred.err.xt,
              timeSpent.String=timeSpent.String))
}