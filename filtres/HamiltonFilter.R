# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Hamilton Filtering
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
# HAMILTON FILTER Implementation du filtre d'Hamilton
#    Permet le calcul des probabilités filtrées : omega_t = P[C_t | X_1:t]
# --------------------------------------------------------
normal.HMM.HamiltonFilter = function(mu,
                                     sigma,
                                     matGamma,
                                     initial.Distribution=NULL,
                                     DATA,
                                     exp.var,
                                     distribution="Normal",
                                     Transition.Type,
                                     type,
                                     nbStepsBack,
                                     nu=4){
  
  nbRegime <- length(mu)
  n <- length(DATA)
  
  
  # -------------------------------------------------
  # HAMILTON FORWARD FILTERING
  # -------------------------------------------------
  u.t <- p.ct.x1t <- p.ct.x1tm1 <- matrix(nrow=nbRegime,ncol=n)
  pred.err.xt <- standard.pred.err.xt <- esp.xt <- std.xt <- matrix(nrow=1,ncol=n)
  diag.Gamma <- matrix(nrow=nbRegime,ncol=n)
  
  if (Transition.Type!="GAS"){
    # Verification that in the event of non-homogeneous transitions, the initial distribution
    #  is correctly specified
    if (!Transition.Type %in% c("Random","Homogeneous")){
      if (is.null(initial.Distribution)){
        stop(paste("HAMILTON FILTER : The initial distribution must be provided explicitly if the Transitions are non-Homogeneous.\nPlease specify a vector of length ",
                   nbRegime," for the 'initial.Distribution' parameter."))
      }
    }
    
    if (is.null(initial.Distribution)){
      if (!type %in% c("UNIVARIEE","MELANGE")){
        initial.Distribution <- solve(t(diag(nbRegime)-matGamma+1),rep(1,nbRegime))
      } else {
        initial.Distribution <- matGamma
      }
    }
    
    if (type %in% c("UNIVARIEE")){
      initial.Distribution <- 1
    }
    
    if (!type %in% c("UNIVARIEE")){
      
      # FIRST STEP OF THE PROCESS
      p.ct.x1tm1[,1] <- initial.Distribution
      u.t[,1] <- initial.Distribution
      esp.xt[1] <- t(u.t[,1]) %*% mu
      std.xt[1] <- sum(t(u.t[,1]) %*% (diag(mu)^2 + diag(sigma)^2))
      
      pred.err.xt[1] <- esp.xt[1] - DATA[1]
      standard.pred.err.xt[1] <- pred.err.xt[1]/std.xt[1]
      
      if (distribution=="Normal"){
        a.j <- p.ct.x1tm1[,1]*(1/sqrt(2*pi*sigma^2))*exp((-(DATA[1]-mu)^2)/(2*sigma^2))
        
      } else if  (distribution=="Student"){
        x <- (DATA[1]-mu)/sigma
        a.j <- p.ct.x1tm1[,1]*((gamma((nu+1)/2)/gamma(nu/2))/sqrt(nu*pi)/((1+(x^2)/nu)^((nu+1)/2)))/sigma
        
      }else {
        stop("The distribution provided is not accepted. \n  Please provide one of the two accepted choices provided here : \n     Normal, Student")
      }
      
      a <- sum(a.j)
      llk <- log(a)
      omega.t <- a.j/a
      p.ct.x1t[,1] <- omega.t
      
      if (!Transition.Type %in% c("LOGIT","LOGIT.w.filter")){
        Gamma.Built <- matGamma
        diag.Gamma <- diag(Gamma.Built)
      }
      
      # On ajoute déjà le dernier pas
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
      
      # STEPS (2:n) OF THE PROCESS
      for (i in 2:n){
        if (Transition.Type %in% c("LOGIT","LOGIT.w.filter")){
          if (Transition.Type=="LOGIT"){
            exp.var.Gamma <- if (is.null(exp.var)){
              DATA[i-1]
            } else {
              exp.var[i-1]
            }
          } else if (Transition.Type=="LOGIT.w.filter"){
            exp.var.Gamma <-if (is.null(exp.var)){
              c(DATA[i-1],p.ct.x1t[,(i-1)])
            } else {
              c(exp.var[i-1],p.ct.x1t[,(i-1)])
            }
          }
          Gamma.Built.obj <- Gamma.Build(prob.i = matGamma,
                                         Transition.Type=Transition.Type,
                                         nbStepsBack=nbStepsBack,
                                         exp.var=exp.var.Gamma,
                                         type=type)
          
          diag.Gamma[,i-1] <- Gamma.Built.obj$diag.Gamma
          Gamma.Built <- Gamma.Built.obj$matGamma
          # print(list(prob.i = matGamma,
          #            Transition.Type=Transition.Type,
          #            nbStepsBack=nbStepsBack,
          #            exp.var=exp.var.Gamma,
          #            type=type))
        }
        
        if (!type %in% c("UNIVARIEE","MELANGE")){
          p.ct.x1tm1[,i] <- omega.t%*%Gamma.Built
          u.t[,i] <- u.t[,i-1] %*% Gamma.Built
        } else if (!type %in% c("UNIVARIEE")){
          p.ct.x1tm1[,i] <- matGamma
          u.t[,i] <- u.t[,i-1]
        } else {
          p.ct.x1tm1[,i] <- 1
          u.t[,i] <- u.t[,i-1]
        }
        
        if (distribution=="Normal"){
          a.j <- p.ct.x1tm1[,i]*((1/sqrt(2*pi*sigma^2))*exp((-(DATA[i]-mu)^2)/(2*sigma^2)))
          
        } else if  (distribution=="Student"){
          x <- (DATA[i]-mu)/sigma
          a.j <- p.ct.x1tm1[,i]*((gamma((nu+1)/2)/gamma(nu/2))/sqrt(nu*pi)/((1+(x^2)/nu)^((nu+1)/2)))/sigma
          
        }
        
        a <- sum(a.j)
        llk <- llk + log(a)
        omega.t <- a.j/a
        p.ct.x1t[,i] <- omega.t
        
        esp.xt[i] <- t(u.t[,i]) %*% mu
        std.xt[i] <- sum(t(u.t[,i]) %*% (diag(mu)^2 + diag(sigma)^2))
        pred.err.xt[i] <- esp.xt[i] - DATA[i]
        standard.pred.err.xt[i] <- pred.err.xt[i]/std.xt[i]
      }
      
    } else {
      llk <- sum(dnorm(DATA, 
                       mean = mu, 
                       sd = sigma, 
                       log = TRUE))
      p.ct.x1tm1 <- matrix(rep(initial.Distribution,n),nrow=1)
      p.ct.x1t  <- matrix(rep(initial.Distribution,n),nrow=1)
      u.t        <- matrix(rep(initial.Distribution,n),nrow=1)
      esp.xt     <- matrix(rep(mu,n),nrow=1)
      std.xt     <- matrix(rep(sum(initial.Distribution * (mu^2 + sigma^2)),n),nrow=1)
      pred.err.xt <- matrix(esp.xt - matrix(DATA,nrow=1),nrow=1)
      standard.pred.err.xt <- pred.err.xt/std.xt
    }
    
    
    return(list(llk=llk,
                p.ct.x1tm1=p.ct.x1tm1,
                p.ct.x1t=p.ct.x1t,
                u.t=u.t,
                esp.xt=esp.xt,
                pred.err.xt=pred.err.xt,
                standard.pred.err.xt=standard.pred.err.xt,
                diag.Gamma=diag.Gamma))
    
  } else if (Transition.Type=="GAS") {
    
    len<-length(DATA)
    GQ<-gauss.quad.prob(30,"normal", mu=median(DATA), sigma=sd(DATA))
    weights<-GQ$weights
    nodes<-GQ$nodes
    
    # print(nodes)
    
    temp<-.C("Filtering_2RegimesGAS",
             "T"=as.integer(len),
             # "initial.Distribution" = as.numeric(initial.Distribution),
             "B"=0,
             "par"=as.numeric(c(mu,
                                sigma,
                                matGamma)), 
             "y"=as.double(DATA),
             "logLikSum"=double(1),
             "X_t"=double(2*len),"X_tlag"=double(2*len),
             "p11"=double(len),"p22"=double(len), 
             "nodes"=as.double(nodes), 
             "weights"=as.double(weights),
             "mu_weights"=mu,
             "sigma_weights"=sigma,
             NAOK=TRUE,PACKAGE="HMMGAS")
    
    
    filter.res <- temp$logLikSum
    
    attr(filter.res,"X.t")    <- matrix(temp$X_t,nrow=2,byrow=F)
    attr(filter.res,"X.tlag") <- matrix(temp$X_tlag,nrow=2,byrow=F)
    attr(filter.res,"p11")<-temp$p11
    attr(filter.res,"p22")<-temp$p22
    u.t <- matrix(temp$X_tlag,nrow=2,byrow=F)
    
    rbind(temp$p1, temp$p2)
    initial.Distribution <- u.t[,1]
    
    esp.xt <- mu %*% u.t
    std.xt <- (mu^2 + sigma^2) %*% u.t
    
    pred.err.xt <- esp.xt - DATA
    standard.pred.err.xt <- pred.err.xt/std.xt
    
    
    return(list(llk=filter.res[1:1],
                p.ct.x1tm1=attr(filter.res, "X.tlag"),
                p.ct.x1t=attr(filter.res, "X.t"),
                initial.Distribution=initial.Distribution,
                u.t=u.t,
                esp.xt=esp.xt,
                pred.err.xt=pred.err.xt,
                standard.pred.err.xt=standard.pred.err.xt,
                diag.Gamma=matrix(c(attr(filter.res, "p11"), 
                                    attr(filter.res, "p22")),
                                  nrow=2, 
                                  byrow=T)))
  }
  
}

