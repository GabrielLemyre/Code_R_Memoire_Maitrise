# ————————————————————————————————————————————————————————————————————————————————————————————————————
# Hidden Markov Models
# PARAMÉTRISATION À ENTRAINER
# ————————————————————————————————————————————————————————————————————————————————————————————————————
# written
# Gabriel LEMYRE
# ————————————————————————————————————————————————————————————————————————————————————————————————————
# Under the supervision of :
# Maciej AUGUSTYNIAK
# ————————————————————————————————————————————————————————————————————————————————————————————————————
# First version : 21 décembre, 2020
# Last version : 21 décembre, 2020
# ————————————————————————————————————————————————————————————————————————————————————————————————————

# ————————————————————————————————————————————————————————
# ENREGISTRE TOUS LES MODÈLES À ENTRAINER DANS UNE LISTE
# ————————————————————————————————————————————————————————
# Retourne l'objet de modèle général
# ————————————————————————————————————————————————————————


# ————————————————————————————————————————————————————————————————————————————————————
# ////////////////////////////////////////////////////////////////////
# SPÉCIFICATIONS À ENTRAÎNER ----
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
# ————————————————————————————————————————————————————————————————————————————————————
Modele.general <- function(variable.explicative,
                           parametres.initiaux){
  parametres.initiaux %>% with({
    return(
      list(
        # ———————————————————————————————————————————————
        # ————————————————————————————————————————————————————————————————————————————————————
        # ////////////////////////////////////////////////////////////////////
        # MODÈLES DE BASES ----
        # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        # ————————————————————————————————————————————————————————————————————————————————————
        modele.Normale = list(type="UNIVARIEE",
                              mu0=mu.Test.Normale,
                              sigma0=sigma.Test.Normale),
        # ———————————————————————————————————————————————
        # ———————————————————————————————————————————————
        modele.Melange.2 = list(type="MELANGE",
                                mu0=mu.Test.2,
                                sigma0=sigma.Test.2,
                                Gamma0=Gamma.Melange.2),
        # ———————————————————————————————————————————————
        modele.Melange.3 = list(type="MELANGE",
                                mu0=mu.Test.3,
                                sigma0=sigma.Test.3,
                                Gamma0=Gamma.Melange.3),
        # ———————————————————————————————————————————————
        modele.Melange.4 = list(type="MELANGE",
                                mu0=mu.Test.4,
                                sigma0=sigma.Test.4,
                                Gamma0=Gamma.Melange.4),
        # ———————————————————————————————————————————————
        # ————————————————————————————————————————————————————————————————————————————————————
        # ////////////////////////////////////////////////////////////////////
        # STRAIGHT FORWARD HMM WITH K=2,3,4 ----
        # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        # ————————————————————————————————————————————————————————————————————————————————————
        modele.HMM.2.10 = list(type="HMM",
                               mu0=mu.Test.2,
                               sigma0=sigma.Test.2,
                               Gamma0=Gamma.HMM.2,
                               initial.Distribution=c(1,0)),
        modele.HMM.2.01 = list(type="HMM",
                               mu0=mu.Test.2,
                               sigma0=sigma.Test.2,
                               Gamma0=Gamma.HMM.2,
                               initial.Distribution=c(0,1)),
        # ———————————————————————————————————————————————
        modele.HMM.3 = list(type="HMM",
                            mu0=mu.Test.3,
                            sigma0=sigma.Test.3,
                            Gamma0=Gamma.HMM.3),
        # ———————————————————————————————————————————————
        # modele.HMM.4.eq = list(type="HMM",
        #                        mu0=mu.Test.4,
        #                        sigma0=sigma.Test.4,
        #                        Gamma0=Gamma.HMM.4.1,
        #                        initial.Distribution=rep(1,4)/4),
        # 
        modele.HMM.4.1000 = list(type="HMM",
                                 mu0=mu.Test.4,
                                 sigma0=sigma.Test.4,
                                 Gamma0=Gamma.HMM.4.1,
                                 initial.Distribution=c(1,0,0,0)),
        # 
        modele.HMM.4.0100 = list(type="HMM",
                                 mu0=mu.Test.4,
                                 sigma0=sigma.Test.4,
                                 Gamma0=Gamma.HMM.4.1,
                                 initial.Distribution=c(0,1,0,0)),
        #
        modele.HMM.4.0010 = list(type="HMM",
                                 mu0=mu.Test.4,
                                 sigma0=sigma.Test.4,
                                 Gamma0=Gamma.HMM.4.1,
                                 initial.Distribution=c(0,0,1,0)),
        #
        modele.HMM.4.0001 = list(type="HMM",
                                 mu0=mu.Test.4,
                                 sigma0=sigma.Test.4,
                                 Gamma0=Gamma.HMM.4.1,
                                 initial.Distribution=c(0,0,0,1)),
        # ————————————————————————————————————————————————————————————————————————————————————
        # ////////////////////////////////////////////////////////////////////
        # HHMM models ----
        # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        # ————————————————————————————————————————————————————————————————————————————————————
        # modele.HHMM.eq = list(type="HHMM",
        #                       mu0=mu.Test.4,
        #                       sigma0=sigma.Test.4,
        #                       Gamma0=Gamma.HHMM.1,
        #                       initial.Distribution=rep(1,4)/4),
        #
        modele.HHMM.1000 = list(type="HHMM",
                                mu0=mu.Test.4,
                                sigma0=sigma.Test.4,
                                Gamma0=Gamma.HHMM.1,
                                initial.Distribution=c(1,0,0,0)),
        #
        modele.HHMM.0100 = list(type="HHMM",
                                mu0=mu.Test.4,
                                sigma0=sigma.Test.4,
                                Gamma0=Gamma.HHMM.2,
                                initial.Distribution=c(0,1,0,0)),
        #
        modele.HHMM.0010 = list(type="HHMM",
                                mu0=mu.Test.4,
                                sigma0=sigma.Test.4,
                                Gamma0=Gamma.HHMM.3,
                                initial.Distribution=c(0,0,1,0)),
        #
        modele.HHMM.0001 = list(type="HHMM",
                                mu0=mu.Test.4,
                                sigma0=sigma.Test.4,
                                Gamma0=Gamma.HHMM.4,
                                initial.Distribution=c(0,0,0,1)),
        # ————————————————————————————————————————————————————————————————————————————————————
        # ////////////////////////////////////////////////////////////////////
        # HMM - DIEBOLD ----
        # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        # ————————————————————————————————————————————————————————————————————————————————————
        # ————————————————————————————————————————————————————————————————————————————————————
        #   VALEURS DÉCALLÉES DE L'INDICE ----
        # ————————————————————————————————————————————————————————————————————————————————————
        # modele.Diebold.eq.past = list(type="HMM",
        #                               Transition.Type="Diebold",
        #                               mu0=mu.Test.4,
        #                               sigma0=sigma.Test.4,
        #                               Gamma0=Gamma.HMM.Diebold.1,
        #                               initial.Distribution=rep(1,4)/4),
        #
        modele.Diebold.1000.past = list(type="HMM",
                                        Transition.Type="Diebold",
                                        mu0=mu.Test.4,
                                        sigma0=sigma.Test.4,
                                        Gamma0=Gamma.HMM.Diebold.1,
                                        initial.Distribution=c(1,0,0,0)),
        #
        modele.Diebold.0100.past = list(type="HMM",
                                        Transition.Type="Diebold",
                                        mu0=mu.Test.4,
                                        sigma0=sigma.Test.4,
                                        Gamma0=Gamma.HMM.Diebold.2,
                                        initial.Distribution=c(0,1,0,0)),
        #
        modele.Diebold.0010.past = list(type="HMM",
                                        Transition.Type="Diebold",
                                        mu0=mu.Test.4,
                                        sigma0=sigma.Test.4,
                                        Gamma0=Gamma.HMM.Diebold.3,
                                        initial.Distribution=c(0,0,1,0)),
        #
        modele.Diebold.0001.past = list(type="HMM",
                                        Transition.Type="Diebold",
                                        mu0=mu.Test.4,
                                        sigma0=sigma.Test.4,
                                        Gamma0=Gamma.HMM.Diebold.4,
                                        initial.Distribution=c(0,0,0,1)),
        # ————————————————————————————————————————————————————————————————————————————————————
        #   VARIABLE EXPLICATIVE SPECIFIÉ EN INTRO ----
        # ————————————————————————————————————————————————————————————————————————————————————
        # modele.Diebold.eq.exp.var = list(type="HMM",
        #                                  Transition.Type="Diebold",
        #                                  exp.var.name=variable.explicative,
        #                                  mu0=mu.Test.4,
        #                                  sigma0=sigma.Test.4,
        #                                  Gamma0=Gamma.HMM.Diebold.1,
        #                                  initial.Distribution=rep(1,4)/4),
        #
        modele.Diebold.1000.exp.var = list(type="HMM",
                                           Transition.Type="Diebold",
                                           exp.var.name=variable.explicative,
                                           mu0=mu.Test.4,
                                           sigma0=sigma.Test.4,
                                           Gamma0=Gamma.HMM.Diebold.1,
                                           initial.Distribution=c(1,0,0,0)),
        #
        modele.Diebold.0100.exp.var = list(type="HMM",
                                           Transition.Type="Diebold",
                                           exp.var.name=variable.explicative,
                                           mu0=mu.Test.4,
                                           sigma0=sigma.Test.4,
                                           Gamma0=Gamma.HMM.Diebold.2,
                                           initial.Distribution=c(0,1,0,0)),
        #
        modele.Diebold.0010.exp.var = list(type="HMM",
                                           Transition.Type="Diebold",
                                           exp.var.name=variable.explicative,
                                           mu0=mu.Test.4,
                                           sigma0=sigma.Test.4,
                                           Gamma0=Gamma.HMM.Diebold.3,
                                           initial.Distribution=c(0,0,1,0)),
        #
        modele.Diebold.0001.exp.var = list(type="HMM",
                                           Transition.Type="Diebold",
                                           exp.var.name=variable.explicative,
                                           mu0=mu.Test.4,
                                           sigma0=sigma.Test.4,
                                           Gamma0=Gamma.HMM.Diebold.4,
                                           initial.Distribution=c(0,0,0,1)),
        # ————————————————————————————————————————————————————————————————————————————————————
        # ////////////////////////////////////////////////////////////////////
        # HMM - GAS ----
        # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        # ————————————————————————————————————————————————————————————————————————————————————
        modele.GAS.1000 = list(type="HMM",
                               Transition.Type="GAS",
                               mu0=mu.Test.2,
                               sigma0=sigma.Test.2,
                               Gamma0=Gamma.HMM.GAS.1),
        #
        modele.GAS.0100 = list(type="HMM",
                               Transition.Type="GAS",
                               mu0=mu.Test.2,
                               sigma0=sigma.Test.2,
                               Gamma0=Gamma.HMM.GAS.2),
        #
        modele.GAS.0010 = list(type="HMM",
                               Transition.Type="GAS",
                               mu0=mu.Test.2,
                               sigma0=sigma.Test.2,
                               Gamma0=Gamma.HMM.GAS.3),
        
        modele.GAS.0001 = list(type="HMM",
                               Transition.Type="GAS",
                               mu0=mu.Test.2,
                               sigma0=sigma.Test.2,
                               Gamma0=Gamma.HMM.GAS.4),
        # ————————————————————————————————————————————————————————————————————————————————————
        # ////////////////////////////////////////////////////////////////////
        # DDMS models ----
        # \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
        # ————————————————————————————————————————————————————————————————————————————————————
        # ————————————————————————————————————————————————————————————————————————————————————
        #   TAU = 3 ----
        # ————————————————————————————————————————————————————————————————————————————————————
        # modele.DDMS.3.eq = list(type="DDMS",
        #                         mu0=mu.Test.4,
        #                         sigma0=sigma.Test.4,
        #                         Gamma0=Gamma.DDMS.3.1,
        #                         nbStepsBack=3,
        #                         initial.Distribution=rep(1,4)/4),
        #
        modele.DDMS.3.1000 = list(type="DDMS",
                                  mu0=mu.Test.4,
                                  sigma0=sigma.Test.4,
                                  Gamma0=Gamma.DDMS.3.1,
                                  nbStepsBack=3,
                                  initial.Distribution=c(1,0,0,0)),
        #
        modele.DDMS.3.0100 = list(type="DDMS",
                                  mu0=mu.Test.4,
                                  sigma0=sigma.Test.4,
                                  Gamma0=Gamma.DDMS.3.2,
                                  nbStepsBack=3,
                                  initial.Distribution=c(0,1,0,0)),
        #
        modele.DDMS.3.0010 = list(type="DDMS",
                                  mu0=mu.Test.4,
                                  sigma0=sigma.Test.4,
                                  Gamma0=Gamma.DDMS.3.3,
                                  nbStepsBack=3,
                                  initial.Distribution=c(0,0,1,0)),
        
        modele.DDMS.3.0001 = list(type="DDMS",
                                  mu0=mu.Test.4,
                                  sigma0=sigma.Test.4,
                                  Gamma0=Gamma.DDMS.3.4,
                                  nbStepsBack=3,
                                  initial.Distribution=c(0,0,0,1)),
        # ————————————————————————————————————————————————————————————————————————————————————
        #   TAU = 5 ----
        # ————————————————————————————————————————————————————————————————————————————————————
        # modele.DDMS.5.eq = list(type="DDMS",
        #                         mu0=mu.Test.4,
        #                         sigma0=sigma.Test.4,
        #                         Gamma0=Gamma.DDMS.5.1,
        #                         nbStepsBack=5,
        #                         initial.Distribution=rep(1,4)/4),
        # 
        modele.DDMS.5.1000 = list(type="DDMS",
                                  mu0=mu.Test.4,
                                  sigma0=sigma.Test.4,
                                  Gamma0=Gamma.DDMS.5.1,
                                  nbStepsBack=5,
                                  initial.Distribution=c(1,0,0,0)),
        # 
        modele.DDMS.5.0100 = list(type="DDMS",
                                  mu0=mu.Test.4,
                                  sigma0=sigma.Test.4,
                                  Gamma0=Gamma.DDMS.5.2,
                                  nbStepsBack=5,
                                  initial.Distribution=c(0,1,0,0)),
        #
        modele.DDMS.5.0010 = list(type="DDMS",
                                  mu0=mu.Test.4,
                                  sigma0=sigma.Test.4,
                                  Gamma0=Gamma.DDMS.5.3,
                                  nbStepsBack=5,
                                  initial.Distribution=c(0,0,1,0)),
        #
        modele.DDMS.5.0001 = list(type="DDMS",
                                  mu0=mu.Test.4,
                                  sigma0=sigma.Test.4,
                                  Gamma0=Gamma.DDMS.5.4,
                                  nbStepsBack=5,
                                  initial.Distribution=c(0,0,0,1)),
        # ————————————————————————————————————————————————————————————————————————————————————
        #   TAU = 20 ----
        # ————————————————————————————————————————————————————————————————————————————————————
        # modele.DDMS.20.eq = list(type="DDMS",
        #                          mu0=mu.Test.4,
        #                          sigma0=sigma.Test.4,
        #                          Gamma0=Gamma.DDMS.20.1,
        #                          nbStepsBack=20,
        #                          initial.Distribution=rep(1,4)/4),
        # 
        modele.DDMS.20.1000 = list(type="DDMS",
                                   mu0=mu.Test.4,
                                   sigma0=sigma.Test.4,
                                   Gamma0=Gamma.DDMS.20.1,
                                   nbStepsBack=20,
                                   initial.Distribution=c(1,0,0,0)),
        # 
        modele.DDMS.20.0100 = list(type="DDMS",
                                   mu0=mu.Test.4,
                                   sigma0=sigma.Test.4,
                                   Gamma0=Gamma.DDMS.20.2,
                                   nbStepsBack=20,
                                   initial.Distribution=c(0,1,0,0)),
        # 
        modele.DDMS.20.0010 = list(type="DDMS",
                                   mu0=mu.Test.4,
                                   sigma0=sigma.Test.4,
                                   Gamma0=Gamma.DDMS.20.3,
                                   nbStepsBack=20,
                                   initial.Distribution=c(0,0,1,0)),
        #
        modele.DDMS.20.0001 = list(type="DDMS",
                                   mu0=mu.Test.4,
                                   sigma0=sigma.Test.4,
                                   Gamma0=Gamma.DDMS.20.4,
                                   nbStepsBack=20,
                                   initial.Distribution=c(0,0,0,1))
      )
    )
  })
}