# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# INFORMATION POUR PARAMÉTRISATION DU S&P TSX
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
# CADUSD - paramétrisation
# --------------------------------------------------------
# 
# --------------------------------------------------------



# TEST SI DOCUMENT APPELÉ DE L'EXTÉRIEUR OU PAS,
#   Si appelé de l'extérieur, 
if(tryCatch(outside.run, error = function(err){return(TRUE)})){
      # Closing all connections
      closeAllConnections()
      # Vide la cache de R
      rm(list=ls())
      # Retrait de tous les graphiques ouvert
      graphics.off()
      outside.run <- TRUE
      ignored.files.OUTSIDE.RUN <- "Parametrisations"
} else {
      ignored.files.OUTSIDE.RUN <- c()
}

# Chemin vers le document qui roule le modèle
path.MAIN <- '~/Documents/GitHub/Memoire_Maitrise/code/r/__modelisation/main.R'

INDPRO=list(
      data.pars = list(index="INDPRO.log",
                       variables.explicatives = c("INDPRO.log"),
                       transformations = c("scale"),
                       scalingFactor = c(100),
                       date.boundaries = c(as.Date("1919-01-01"),as.Date("2020-12-26")), # Partagé par toutes les modèles
                       freq.str="M"
      ),
      general.pars = list(
            # type.a.rouler = c("MELANGE","HMM"),
            load.model = "/Users/gabriellemyre/Documents/GitHub/Memoire_Maitrise/__resultats/__Résumé résultats pour Maciej/INDPRO/[ 1919-01-01 ~ 2020-12-26 ]/INDPROlog - (M) - [1919-01-01, 2020-12-26] - 23h37m18s.rds",
            # test.Transformations=TRUE,
            # 
            mu0 = list(
                  "1" = c(0.03),
                  "2" = c(0.19, 0.17),
                  "3" = c(0.16,  0.3,  -0.0513),
                  "4" = c(-0.606,  1.8538,  0.2513,  0.4)
            ),
            # 
            sigma0 = list(
                  "1" = c(1.607599),
                  "2" = c(2, 1.5),
                  "3" = c(3, 2, 1),
                  "4" = c(4, 3, 2, 1)
            )
      ),
      UNIVARIEE = list(),
      # 
      MELANGE = list(
            K = c(2,
                  3,
                  4),
            gamma0 = list("2" = c(0.64, 0.36),
                          "3" = c(0.20, 0.44, 0.36),
                          "4" = c(0.19, 0.11, 0.43, 0.36))
      ),
      # 
      HMM = list(
            K = c(2,
                  3,
                  4),
            Transition.Type = "Homogeneous",
            gamma0 = list(
                  "2" = matrix(c(0.981103103, 0.0188969, 
                                 0.008198487, 0.9918015),
                               ncol=2 , byrow=T),
                  "3" = matrix(c(0.98, 0.01, 0.01,
                                 0.01, 0.98, 0.01,
                                 0.01, 0.01, 0.98),
                               byrow=T,nrow=3),
                  "4" = matrix(c(0.97, 0.01, 0.01, 0.01, 
                                 0.01, 0.97, 0.01, 0.01,
                                 0.01, 0.01, 0.97, 0.01,
                                 0.01, 0.01, 0.01, 0.97),
                               ncol=4 , byrow=T)
            ),
            initial.Distribution = list(
                  "2" = list(c(0.25, 0.75)),
                  "3" = list(c(0,1,0)),
                  "4" = list(c(0, 1, 0, 0)))
      ),
      # 
      HHMM = list(
            K=4,
            gamma0 = list(
                  "4" = matrix(
                        c(0.97, 0.02, 0.01,
                          0.97, 0.02, 0.01,
                          0.97, 0.02, 0.01,
                          0.97, 0.02, 0.01,
                          0.2,       0.8,          0,
                          0.9,       0.1,          0,
                          0.4,       0.6,          0,
                          0.6,       0.4,          0),
                        ncol=3,byrow=T)
            ),
            initial.Distribution = list(
                  "4" = list(
                        c(0.051, 0.57, 0.373, 0.006)
                  )
            )
      ),
      # 
      LOGIT = list(
            K=4,
            Transition.Type = "LOGIT",
            gamma0 = list(
                  "SP.500" = c(11.9153937,  -3.2901744,  -3.2196457,
                               9.6929927,   4.6357821,   3.6808849,
                               4.8872555,   9.7381110,  -2.8478310,
                               -11.7162761,   5.8048835,   1.0552933,
                               3.3021530,   1.9993717,   1.3619955,
                               -0.4021961,   0.2381303,  -0.4848302,
                               -1.5521677,  -0.5938818,  -3.0357197,
                               -1.9905202,  -3.1908825,  -7.6027283),
                  "INDPRO" = c(12.0525027,  -3.2897033,  -3.2194887,
                                   25.3769695,  21.8681980,  20.1045678,
                                   21.8813330,  25.8046022, -10.2094039,
                                   -13.0967651,  18.2532958,  21.7022517,
                                   1.6461908,   1.9992837,   1.3620767,
                                   0.0836822,   1.4733259,  -0.9051992,
                                   4.9144567,   3.6489609,  -8.7188818,
                                   4.1613465,  -0.6103364,  -2.4674836),
                  "CADUSD" = c(11.9153937,  -3.2901744,  -3.2196457,
                               9.6929927,   4.6357821,   3.6808849,
                               4.8872555,   9.7381110,  -2.8478310,
                               -11.7162761,   5.8048835,   1.0552933,
                               3.3021530,   1.9993717,   1.3619955,
                               -0.4021961,   0.2381303,  -0.4848302,
                               -1.5521677,  -0.5938818,  -3.0357197,
                               -1.9905202,  -3.1908825,  -7.6027283)
            ),
            initial.Distribution = list(
                  "4" = list(
                        # c(0, 1, 0, 0)
                        c(0.051, 0.57, 0.373, 0.006)
                  )
            )
      ),
      #
      GAS = list(
            K=2,
            Transition.Type = "GAS",
            gamma0 = list(
                  "2" = c(0.25, 0.75, # Distribution initiale
                          0.35, 0.99, # Diagonale de la matrice A
                          0.7, 0.7)
            )
      ),
      # 
      DDMS = list(
            K = 4,
            tau = c(
                  3,
                  5,
                  14
            ),
            gamma0 = list(
                  "3" = c(
                    8.96, -1.03, -2.92,
                    4.09,  4.85,  2.05,
                    4.96,  5.86,  3.87,
                    -4.52,  4.13, 1.82,
                    1.69, -0.80,  2.08,
                    0.54,  0.26, -0.30,
                    6.38,  5.789, -0.1,
                    -1.01, -2.07,  1.5
                  ),
                  "5" = c(
                    8.96, -1.03, -2.92,
                    4.09,  4.85,  2.05,
                    4.96,  5.86,  3.87,
                    -4.52,  4.13, 1.82,
                    1.69, -0.80,  2.08,
                    0.54,  0.26, -0.30,
                    6.38,  5.789, -0.1,
                    -1.01, -2.07,  1.5
                  ),
                  "14" = c(
                    8.96, -1.03, -2.92,
                    4.09,  4.85,  2.05,
                    4.96,  5.86,  3.87,
                    -4.52,  4.13, 1.82,
                    1.69, -0.80,  2.08,
                    0.54,  0.26, -0.30,
                    6.38,  5.789, -0.1,
                    -1.01, -2.07,  1.5
                  )
            ),
            initial.Distribution = list(
                  "4" = list(
                        c(0, 1, 0, 0)
                  )
            )
      )
)

if(outside.run){
      cat("Running", path.MAIN)
      run.model <- INDPRO
      source(path.MAIN)
}


