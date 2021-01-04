# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# TEST - Unconstraining and reconstraining parameters
# ----------------------------------------------------------------------------------------------------
# written
# Gabriel LEMYRE
# ----------------------------------------------------------------------------------------------------
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ----------------------------------------------------------------------------------------------------
# Last version : 22 avril 2020
# Last version : 22 avril 2020
# ----------------------------------------------------------------------------------------------------

# --------------------------------------------------------
# PARAMÈTRES CONTRAINTS -> PARAMÈTRES NON-CONTRAINTS -> PARAMÈTRES CONTRAINTS
# --------------------------------------------------------
# Function to test the conversion from natural parameters 
# to their working counterparts
#     The idea is that by changing the constrained parameters to
#     unconstrained versions, the optimization can be done without
#     constraints
# --------------------------------------------------------
N2N.test = function(mu, sigma, 
                    matGamma=NULL,
                    type, 
                    Transition.Type=NULL,
                    print.steps=FALSE){
  
  input <- list(mu=mu,
                sigma=sigma,
                matGamma=matGamma)
  
  if (print.steps){
    hline()
    cat("Input :\n")
    hline()
    print(input)
  }
  
  # retrait des contraintes sur les paramètres
  parvect.test <- normal.HMM.N2W(mu,sigma,matGamma,type=type,Transition.Type=Transition.Type)
  if (print.steps){
    hline()
    cat("parvect.test :\n")
    hline()
    print(parvect.test)
  }
  
  
  # retour aux version contraintes des paramètres
  natpar.test <- normal.HMM.W2N(parvect.test,type=type,Transition.Type=Transition.Type)
  if (print.steps){
    hline()
    cat("natpar.test :\n")
    hline()
    print(natpar.test)
  }
  
  # Message de début du test
  hline()
  cat("Erreurs maximales lors de la conversion de contraints -> non-contraints -> contraints\n")
  hline() # Impression ligne standard de séparation dans la console
  
  if (type!="UNIVARIEE"){
    err <- diff.liste(input,natpar.test)
  } else {
    err <- diff.liste(input[!names(input) %in% c("matGamma")],
                      natpar.test[!names(natpar.test) %in% c("matGamma")])
  }
  
  
  hline() # Impression ligne standard de séparation dans la console
  
  # retourne le paramètres non-contraint
  return(err)
}
  
  
  
  
  
  
  