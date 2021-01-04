# ----------------------------------------------------------------------------------------------------
# Hidden Markov Models
# INFORMATION POUR PARAMÉTRISATION DU S&P 500
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

# # Closing all connections
# closeAllConnections()
# # Vide la cache de R
# rm(list=ls())
# # Retrait de tous les graphiques ouvert
# graphics.off()
# --------------------------------------------------------
# S&P 500 - paramétrisation
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

SP.500=list(
  data.pars = list(index="SP.500",
                   variables.explicatives = c(
                     "INDPRO.log",
                     "SP.500"
                   ),
                   transformations = c(
                     "logR",
                     "scale"
                   ),
                   scalingFactor = c(
                     100,
                     100
                   ),
                   date.boundaries = c(as.Date("1999-12-31"),Sys.Date()), # Partagé par toutes les modèles
                   freq.str="Q"
  ),
  general.pars = list(
    type.a.rouler = c(
      "HMM",
      "DDMS",
      "LOGIT"
    ),
    # load.model = "/Users/gabriellemyre/Documents/GitHub/Memoire_Maitrise/__resultats/__FINAUX__/SP500/[ 1999-12-31 ~ 2020-12-26 ]/SP500 - (Q) - [1999-12-31, 2020-12-26] - FINAUX.rds",
    # test.Transformations=TRUE,
    # 
    mu0 = list(
      "1" = c(0.03),
      "2" = c(-0.16,  0.06),
      "3" = c(-0.7, 0.2,  0.6),
      "4" = c(-0.27, -0.08,  -0.002,  0.09)
    ),
    # 
    sigma0 = list(
      "1" = c(2.5),
      "2" = c(2.227, 0.731),
      "3" = c(4.5,  1.7, 0.38),
      "4" = c(3.51, 1.56, 0.95, 0.51)
    )
  ),
  UNIVARIEE = list(),
  # 
  MELANGE = list(
    K = c(2,
          3,
          4),
    gamma0 = list("2" = c(0.25, 0.75),
                  "3" = c(0.027, 0.478, 0.495),
                  "4" = c(0.051, 0.570, 0.373, 0.006))
  ),
  # 
  HMM = list(
    K = c(
      # 2,
      # 3,
      4
    ),
    Transition.Type = "Homogeneous",
    gamma0 = list(
      "2" = matrix(c(0.981103103, 0.0188969, 
                     0.008198487, 0.9918015),
                   ncol=2 , byrow=T),
      "3" = matrix(c(0.98, 0.01, 0.01,
                     0.01, 0.98, 0.01,
                     0.01, 0.01, 0.98),
                   byrow=T,nrow=3),
      "4" = matrix(c(0.96827, 0.03171, 0.00001, 0.00001,
                     0.03632, 0.95325, 0.00976, 0.00067,
                     0.00001, 0.01527, 0.97914, 0.00558,
                     0.00001, 0.00001, 0.03702, 0.96296),
                   ncol=4 , byrow=T)
    ),
    initial.Distribution = list(
      "2" = list(c(0.25, 0.75)),
      "3" = list(c(0,1,0)),
      "4" = list(c(0, 1, 0, 0))
    )
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
          0.9,       0.1,          0,
          0.9,       0.1,          0,
          0.6,       0.4,          0,
          0.6,       0.4,          0),
        ncol=3,byrow=T)
    ),
    initial.Distribution = list(
      "4" = list(
        c(0, 1, 0, 0)
      )
    )
  ),
  # 
  HHMM.simplifie = list(
    K=4,
    gamma0 = list(
      "4" = matrix(
        c(0.98, 0.02, # Matrice de tansition parents
          0.02, 0.98,
          0.98, 0.02,  # Matrice de tansition enfants 1
          0.02, 0.98,
          0.98, 0.02,  # Matrice de tansition enfants 2
          0.02, 0.98),
        ncol=2,byrow=T)
    ),
    initial.Distribution = list(
      "4" = list(
        c(0, 1, 0, 0)
      )
    )
  ),
  # 
  LOGIT = list(
    K=4,
    Transition.Type = "LOGIT",
    gamma0 = list(
      # "SP.500" = c( # De INDPRO qui a convergé (Minimise l'erreur de prévision)
      #   5.004786, -118.439484,  -30.083290,
      #   24.191207,   20.692635,   19.105863,
      #   44.211203,   47.625005,  -54.777550,
      #   -13.096765,   18.032581,   21.875960,
      #   0.415372,  -27.057729,   -4.308149,
      #   5.559103,    5.441797,    5.612090,
      #   15.263068,   16.201257,  -30.994257,
      #   4.161346,   -1.774344,   -1.329762
      # ),
      # "SP.500" = c( # De SP.500 qui a convergé (Minimise le BIC)
      #   16.3988563,  12.9266050,  -4.6239321,
      #   6.7280885,   2.1738118,   1.8086030,
      #   12.7942859,  16.7401119,   8.7226543,
      #   -11.7162130,  13.9087757,   1.1149435,
      #   0.4819935,   0.7303530,   0.0604108,
      #   1.7024309,   2.7614406,   2.0216466,
      #   -0.8110554,   0.7464651,  -1.8700124,
      #   -1.9904603,   4.3098921,  -7.6568734
      # ),
      "SP.500" = c( # De SP.500 qui a convergé (Minimise le BIC)
        1,  2,  -2,
        1,   1,   2,
        1,  2,   2,
        -1,  2,   2,
        2,   1,   1,
        1,   1,   2,
        -1,   1,  -2,
        -2,   2,  -2
      ),
      "INDPRO.log" = c( # De SP.500 qui a convergé (Minimise le BIC)
        1,  2,  -2,
        1,   1,   2,
        1,  2,   2,
        -1,  2,   2,
        2,   1,   1,
        1,   1,   2,
        -1,   1,  -2,
        -2,   2,  -2
      ),
      #   "INDPRO.log" = c(12.0525027,  -3.2897033,  -3.2194887,
      #                    25.3769695,  21.8681980,  20.1045678,
      #                    21.8813330,  25.8046022, -10.2094039,
      #                    -13.0967651,  18.2532958,  21.7022517,
      #                    1.6461908,   1.9992837,   1.3620767,
      #                    0.0836822,   1.4733259,  -0.9051992,
      #                    4.9144567,   3.6489609,  -8.7188818,
      #                    4.1613465,  -0.6103364,  -2.4674836),
      #   "CADUSD" = c(11.9153937,  -3.2901744,  -3.2196457,
      #                9.6929927,   4.6357821,   3.6808849,
      #                4.8872555,   9.7381110,  -2.8478310,
      #                -11.7162761,   5.8048835,   1.0552933,
      #                3.3021530,   1.9993717,   1.3619955,
      #                -0.4021961,   0.2381303,  -0.4848302,
      #                -1.5521677,  -0.5938818,  -3.0357197,
      #                -1.9905202,  -3.1908825,  -7.6027283)
      # ),
      initial.Distribution = list(
        "4" = list(
          # c(0.051, 0.570, 0.373, 0.006),
          c(0, 1, 0, 0)
          # c(0, 0, 1, 0)
        )
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
        1,  -2, -2,
        2,   2 ,  1 ,
        -2 , -2 ,-2 ,
        -2,  2,   2,
        -2 , -2 ,  2,
        2,  2,  2 ,
        1 ,  1, -1,
        -1, -2,   3
      ),
      "5" = c(
        2,  -1,  -2,
        2,   2,  -2,
        -2,  -1, -2,
        -2,  1,   1,
        -1,  -2,   1,
        2,   2,   1,
        1,   1, -2,
        -1,  -2,   3
      ),
      "14" = c(
        2,  -1,  -2,
        2,   2,  -2,
        -2,  -1, -2,
        -2,  1,   1,
        -1,  -2,   1,
        2,   2,   1,
        1,   1, -2,
        -1,  -2,   3
      )
    ),
    initial.Distribution = list(
      "4" = list(
        # c(1, 0, 0, 0),
        c(0, 1, 0, 0)
      )
    )
  )
)

if(outside.run){
  cat("Running", path.MAIN)
  run.model <- SP.500
  source(path.MAIN)
}


