# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# Construction de la matrice Gamma
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
#   Construction de la matrice de transition (Gamma)
# --------------------------------------------------------
#   Cette fonction construit la matrice de transition d'un modèle HHMM à 2
#   états parents avec 2 états enfants chaques et d'un modèle FHMM par le 
#   biais des probabilités de transitions horizontales (p1-p2) et transitions
#   verticales (q1-q2)
# -------------------------------
# STRUCTURE POUR LE HIERARCHICAL HMM
# -------------------------------
#     Gamma.i = [a^2.11  +  a^2.12  +  e1     ;  = 1
#                a^2.22  +  a^2.21  +  e2     ;  = 1
#                a^2.33  +  a^2.34  +  e3     ;  = 1
#                a^2.44  +  a^2.43  +  e4     ;  = 1
#                a^1.11  +  a^1.12  +  0      ;  = 1
#                a^1.22  +  a^1.21  +  0      ;  = 1
#                pi.1    +  pi.2    +  0      ;  = 1
#                pi.3    +  pi.4    +  0    ] ;  = 1
# -------------------------------
# STRUCTURE POUR LE FACTORIAL HMM
# -------------------------------
# #          p1          (1-p1)          (1-p2)        p2
# Gamma.1 = [prob.i(1,1) 1-prob.i(1,1) ; 1-prob.i(1,2) prob.i(1,2)];
# #          q1          (1-q1)          (1-q2)        q2
# Gamma.2 = [prob.i(1,3) 1-prob.i(1,3) ; 1-prob.i(1,4) prob.i(1,4)];
# # Produit de Kronecker entre les 2 matrices
# Gamma=kron(Gamma.1,Gamma.2);
# --------------------------------------------------------

Gamma.Build <- function(prob.i, type = NULL,
                       Transition.Type="Homogeneous",
                       nbStepsBack=0, 
                       exp.var=NULL){
  
  matGamma = NULL
  diag.Gamma = NULL
  beta0 <- beta1 <-  NULL
  
  if (Transition.Type=="Homogeneous"){
    if (type=="HHMM"){ # Récupération des paramètres pour HHMM
      a2.11 <- prob.i[1,1]
      a2.12 <- prob.i[1,2]
      e1 <- prob.i[1,3]
      a2.22 <- prob.i[2,1]
      a2.21 <- prob.i[2,2]
      e2 <- prob.i[2,3]
      a2.33 <- prob.i[3,1]
      a2.34 <- prob.i[3,2]
      e3 <- prob.i[3,3]
      a2.44 <- prob.i[4,1]
      a2.43 <- prob.i[4,2]
      e4 <- prob.i[4,3]
      a1.11 <- prob.i[5,1]
      a1.12 <- prob.i[5,2]
      a1.22 <- prob.i[6,1]
      a1.21 <- prob.i[6,2]
      pi.1 <- prob.i[7,1]
      pi.2 <- prob.i[7,2]
      pi.3 <- prob.i[8,1]
      pi.4 <- prob.i[8,2]
      
      # Construction de la matrice de transition à partir de ces paramètres
      Gamma.Built = matrix(c(a2.11+e1*a1.11*pi.1,  a2.12+e1*a1.11*pi.2,          e1*a1.12*pi.3,        e1*a1.12*pi.4,
                             a2.21+e2*a1.11*pi.1,  a2.22+e2*a1.11*pi.2,          e2*a1.12*pi.3,        e2*a1.12*pi.4,
                             e3*a1.21*pi.1,        e3*a1.21*pi.2,          a2.33+e3*a1.22*pi.3,  a2.34+e3*a1.22*pi.4,
                             e4*a1.21*pi.1,        e4*a1.21*pi.2,          a2.43+e4*a1.22*pi.3,  a2.44+e4*a1.22*pi.4),
                           nrow=4,byrow=T)
      
    } else if (type =="HHMM.simplifie"){
      
      Gamma.Parent = prob.i[1:2,1:2]
      Gamma.Enfant.1 = prob.i[3:4,1:2]
      Gamma.Enfant.2 = prob.i[5:6,1:2]
      Gamma.Built = cbind(Gamma.Parent[,1] %x% Gamma.Enfant.1,
                          Gamma.Parent[,2] %x% Gamma.Enfant.2)
      
    } else if (type=="FHMM"){
      #          p1          (1-p1)          (1-p2)        p2
      Gamma.1 = matrix(c(prob.i(1,1), 1-prob.i(1,1), 1-prob.i(1,2), prob.i(1,2)),nrow=2,byrow=T)
      #          q1          (1-q1)          (1-q2)        q2
      Gamma.2 = matrix(c(prob.i(1,3), 1-prob.i(1,3), 1-prob.i(1,4), prob.i(1,4)),nrow=2,byrow=T)
      # Produit de Kronecker entre les 2 matrices
      Gamma.Built = kronecker(Gamma.1,Gamma.2);
      
    } else if (type=="DDMS"){
      # prob.i <- matrix(prob.i,nrow=4,byrow=T)
      
      maxDDMS <- nbStepsBack
      n.coef <- length(prob.i)
      
      
      nbRegime <- (1+sqrt(1+4*(n.coef/2)))/2
      
      beta0 <- matrix(0,nrow=nbRegime,ncol=nbRegime)
      beta1 <- matrix(0,nrow=nbRegime,ncol=nbRegime)
      
      beta0[1,1:(nbRegime-1)] <- prob.i[1:(nbRegime-1)]
      beta0[2:nbRegime,2:nbRegime] <- matrix(prob.i[nbRegime:(n.coef/2)], ncol=nbRegime-1,byrow=T)
      
      prob.i.half <- prob.i[-c(1:(n.coef/2))]
      beta1[1,1:(nbRegime-1)] <- prob.i.half[1:(nbRegime-1)]
      beta1[2:nbRegime,2:nbRegime] <- matrix(prob.i.half[nbRegime:(n.coef/2)], ncol=nbRegime-1,byrow=T)
      
      
      dim.Gamma <- nbRegime*maxDDMS
      matGamma <- matrix(1e-16,ncol=dim.Gamma,nrow=dim.Gamma)
      for (i in 1:dim.Gamma){
        for (j in 1:dim.Gamma){
          istar <- ceiling(i/maxDDMS)
          jstar <- ceiling(j/maxDDMS)
          di <- i-maxDDMS*(istar-1)
          dj <- j-maxDDMS*(jstar-1)
          
          if ((istar==jstar & dj==min(di+1,maxDDMS)) | (istar!=jstar & dj==1)){
            matGamma[i,j] <- exp(beta0[istar,jstar] + beta1[istar,jstar] * di)
          }
        }
      }
      
      row.sum <- as.numeric(matGamma %*% matrix(rep(1,dim.Gamma),ncol=1))
      # row.sum[row.sum==Inf] <- max(row.sum[row.sum!=Inf])
      Gamma.Built <- matGamma / row.sum
    }
    
  } else {
    if (Transition.Type=="LOGIT"){
      
      if (is.null(exp.var)){stop("Please provide the past observation for the non-homogeneous transition probabilities model.")}

      n.coef <- length(prob.i)
      
      nbRegime <- (1+sqrt(1+4*(n.coef/2)))/2
      
      beta0 <- matrix(0,nrow=nbRegime,ncol=nbRegime)
      beta1 <- matrix(0,nrow=nbRegime,ncol=nbRegime)
      
      beta0[1,1:(nbRegime-1)] <- prob.i[1:(nbRegime-1)]
      beta0[2:nbRegime,2:nbRegime] <- matrix(prob.i[nbRegime:(n.coef/2)], ncol=(nbRegime-1),byrow=T)
      
      prob.i.half <- prob.i[-c(1:(n.coef/2))]
      beta1[1,1:(nbRegime-1)] <- prob.i.half[1:(nbRegime-1)]
      beta1[2:nbRegime,2:nbRegime] <- matrix(prob.i.half[nbRegime:(n.coef/2)], ncol=(nbRegime-1),byrow=T)
      
      Gamma.Built.noexp <- beta0 + beta1 * exp.var
      
      Gamma.Built <- exp(Gamma.Built.noexp) / c(exp(Gamma.Built.noexp)%*%rep(1,dim(Gamma.Built.noexp)[2]))
    }
    if (Transition.Type=="LOGIT.final"){
    
      n.coef <- length(prob.i)
      
      nbRegime <- (1+sqrt(1+4*(n.coef/2)))/2
      
      beta0 <- matrix(0,nrow=nbRegime,ncol=nbRegime)
      beta1 <- matrix(0,nrow=nbRegime,ncol=nbRegime)
      
      beta0[1,1:(nbRegime-1)] <- prob.i[1:(nbRegime-1)]
      beta0[2:nbRegime,2:nbRegime] <- matrix(prob.i[nbRegime:(n.coef/2)], ncol=nbRegime-1,byrow=T)
      
      prob.i.half <- prob.i[-c(1:(n.coef/2))]
      beta1[1,1:(nbRegime-1)] <- prob.i.half[1:(nbRegime-1)]
      beta1[2:nbRegime,2:nbRegime] <- matrix(prob.i.half[nbRegime:(n.coef/2)], ncol=nbRegime-1,byrow=T)
      
      return(list(beta0=beta0,
                  beta1=beta1))
    }
  }
  
  return(list(matGamma=Gamma.Built,
              diag.Gamma=diag(Gamma.Built),
              beta0=beta0,
              beta1=beta1))
}